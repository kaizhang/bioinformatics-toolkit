{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DeriveGeneric #-}
--------------------------------------------------------------------------------
-- |
-- Module      :  $Header$
-- Copyright   :  (c) 2014 Kai Zhang
-- License     :  MIT

-- Maintainer  :  kai@kzhang.org
-- Stability   :  experimental
-- Portability :  portable

-- functions for processing BED files
--------------------------------------------------------------------------------

module Bio.Data.Bed
    ( BEDFormat(..)
    , BED(..)
    , BED3(..)
    , fetchSeq
    , Sorted(..)
    , compareBed
    , sortBed
    ) where

import Bio.Seq
import Bio.Seq.IO
import Bio.Utils.Misc (readInt, readDouble)
import Control.Monad.State.Strict
import Control.Monad.ST
import qualified Data.ByteString.Char8 as B
import Data.Conduit
import qualified Data.Conduit.List as CL
import Data.Default.Generics
import Data.Maybe
import qualified Data.Vector as V
import qualified Data.Vector.Algorithms.Intro as I
import GHC.Generics
import System.IO

class BEDFormat b where
    fromLine :: B.ByteString -> b
    toLine :: b -> B.ByteString
    -- | field accessor
    chrom :: b -> B.ByteString
    chromStart :: b -> Int
    chromEnd :: b -> Int

    size :: b -> Int
    size bed = chromEnd bed - chromStart bed + 1

    hReadBed :: Handle -> Source IO b
    hReadBed h = do eof <- liftIO $ hIsEOF h
                    unless eof $ do
                        line <- liftIO $ B.hGetLine h
                        yield $ fromLine line
                        hReadBed h
    {-# INLINE hReadBed #-}

    -- | non-streaming version
    hReadBed' :: Handle -> IO [b]
    hReadBed' h = hReadBed h $$ CL.consume
    {-# INLINE hReadBed' #-}

    readBed :: FilePath -> Source IO b
    readBed fl = do handle <- liftIO $ openFile fl ReadMode
                    hReadBed handle
                    liftIO $ hClose handle
    {-# INLINE readBed #-}

    -- | non-streaming version
    readBed' :: FilePath -> IO [b]
    readBed' fl = readBed fl $$ CL.consume
    {-# INLINE readBed' #-}

    writeBed :: FilePath -> [b] -> IO ()
    writeBed fl beds = withFile fl WriteMode $ \h -> mapM_ (B.hPutStrLn h.toLine) beds
    {-# INLINE writeBed #-}

    {-# MINIMAL fromLine, toLine, chrom, chromStart, chromEnd #-}

-- | BED6 format, as described in http://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
data BED = BED
    { _chrom :: !B.ByteString
    , _chromStart :: {-# UNPACK #-} !Int
    , _chromEnd :: {-# UNPACK #-} !Int
    , _name :: !(Maybe B.ByteString)
    , _score :: !(Maybe Double)
    , _strand :: !(Maybe Bool)  -- ^ True: "+", False: "-"
    } deriving (Show, Generic)

instance Default BED

instance BEDFormat BED where
    fromLine l = evalState (f (B.split '\t' l)) 1
      where
        f :: [B.ByteString] -> State Int BED
        f [] = do i <- get
                  if i <= 3 then error "Read BED fail: Incorrect number of fields"
                            else return def
        f (x:xs) = do 
            i <- get
            put (i+1)
            bed <- f xs
            case i of
                1 -> return $ bed {_chrom = x}
                2 -> return $ bed {_chromStart = readInt x}
                3 -> return $ bed {_chromEnd = readInt x}
                4 -> return $ bed {_name = guard' x}
                5 -> return $ bed {_score = getScore x}
                6 -> return $ bed {_strand = getStrand x}
                _ -> return def

        guard' x | x == "." = Nothing
                 | otherwise = Just x
        getScore x | x == "." = Nothing
                   | otherwise = Just . readDouble $ x
        getStrand str | str == "-" = Just False
                      | str == "+" = Just True
                      | otherwise = Nothing
    {-# INLINE fromLine #-}

    toLine (BED f1 f2 f3 f4 f5 f6) = B.intercalate "\t" [ f1
                                                        , (B.pack.show) f2
                                                        , (B.pack.show) f3
                                                        , fromMaybe "." f4
                                                        , score'
                                                        , strand'
                                                        ]
      where
        strand' | f6 == Just True = "+"
                | f6 == Just False = "-"
                | otherwise = "."
        score' = case f5 of
                     Just x -> (B.pack.show) x
                     _ -> "."
    {-# INLINE toLine #-}

    chrom = _chrom
    chromStart = _chromStart
    chromEnd = _chromEnd

data BED3 = BED3 !B.ByteString !Int !Int deriving (Show, Generic)

instance Default BED3 

instance BEDFormat BED3 where
    fromLine l = case B.split '\t' l of
                    (a:b:c:_) -> BED3 a (readInt b) $ readInt c
                    _ -> error "Read BED fail: Incorrect number of fields"
    {-# INLINE fromLine #-}
    
    toLine (BED3 f1 f2 f3) = B.intercalate "\t" [f1, (B.pack.show) f2, (B.pack.show) f3]
    {-# INLINE toLine #-}

    chrom (BED3 f1 _ _) = f1
    chromStart (BED3 _ f2 _) = f2
    chromEnd (BED3 _ _ f3) = f3

fetchSeq :: BioSeq DNA a => Genome -> Conduit BED IO (DNA a)
fetchSeq g = do gH <- liftIO $ gHOpen g
                table <- liftIO $ readIndex gH
                conduitWith gH table
                liftIO $ gHClose gH
  where
    conduitWith h index' = do 
        bed <- await
        case bed of
            Just (BED chr start end _ _ isForward) -> do 
                dna <- liftIO $ getSeq h index' (chr, start, end)
                case isForward of
                    Just False -> yield $ rc dna
                    _ -> yield dna
                conduitWith h index'
            _ -> return ()
{-# INLINE fetchSeq #-}

-- | a type to imply that underlying data is sorted
newtype Sorted b = Sorted b

compareBed :: (BEDFormat b1, BEDFormat b2) => b1 -> b2 -> Ordering
compareBed x y = compare x' y'
  where
    x' = (chrom x, chromStart x, chromEnd x)
    y' = (chrom y, chromStart y, chromEnd y)
{-# INLINE compareBed #-}

-- | sort BED, first by chromosome (alphabetical order), then by chromStart, last by chromEnd
sortBed :: BEDFormat b => V.Vector b -> Sorted (V.Vector b)
sortBed beds = Sorted $ runST $ do
    v <- V.unsafeThaw beds
    I.sortBy compareBed v
    V.unsafeFreeze v
{-# INLINE sortBed #-}
