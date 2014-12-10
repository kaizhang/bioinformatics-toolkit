{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE BangPatterns #-}
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
    ( BEDLike(..)
    , splitBed
    , splitBedBySize
    , Sorted(..)
    , sortBed
    , mergeBed
    , mergeBed'
    -- * BED6 format
    , BED(..)
    -- * BED3 format
    , BED3(..)
    , fetchSeq
    , compareBed
    ) where

import Bio.Seq
import Bio.Seq.IO
import Bio.Utils.Misc
import Control.Monad.State.Strict
import Control.Monad.ST
import qualified Data.ByteString.Char8 as B
import Data.Conduit
import qualified Data.Conduit.List as CL
import Data.Default.Class
import Data.Maybe
import qualified Data.Vector as V
import qualified Data.Vector.Algorithms.Intro as I
import GHC.Generics
import System.IO

-- | a class representing BED-like data, e.g., BED3, BED6 and BED12. BED format
-- uses 0-based index.
class BEDLike b where
    -- | construct bed record from chromsomoe, start location and end location
    asBed :: B.ByteString -> Int -> Int -> b

    -- | convert bytestring to bed format
    fromLine :: B.ByteString -> b

    -- | convert bed to bytestring
    toLine :: b -> B.ByteString
    
    -- | field accessor
    chrom :: b -> B.ByteString
    chromStart :: b -> Int
    chromEnd :: b -> Int

    -- | return the size of a bed region
    size :: b -> Int
    size bed = chromEnd bed - chromStart bed + 1

    hReadBed :: Handle -> Source IO b
    hReadBed h = do eof <- lift $ hIsEOF h
                    unless eof $ do
                        line <- lift $ B.hGetLine h
                        yield $ fromLine line
                        hReadBed h
    {-# INLINE hReadBed #-}

    -- | non-streaming version
    hReadBed' :: Handle -> IO [b]
    hReadBed' h = hReadBed h $$ CL.consume
    {-# INLINE hReadBed' #-}

    readBed :: FilePath -> Source IO b
    readBed fl = do handle <- lift $ openFile fl ReadMode
                    hReadBed handle
                    lift $ hClose handle
    {-# INLINE readBed #-}

    -- | non-streaming version
    readBed' :: FilePath -> IO [b]
    readBed' fl = readBed fl $$ CL.consume
    {-# INLINE readBed' #-}

    writeBed :: FilePath -> Sink b IO ()
    writeBed fl = do handle <- lift $ openFile fl WriteMode
                     go handle
      where
        go h = do x <- await
                  case x of
                      Nothing -> lift $ hClose h
                      Just bed -> (lift . B.hPutStrLn h . toLine) bed >> go h
    {-# INLINE writeBed #-}

    writeBed' :: FilePath -> [b] -> IO ()
    writeBed' fl beds = withFile fl WriteMode $ \h -> mapM_ (B.hPutStrLn h.toLine) beds
    {-# INLINE writeBed' #-}

    {-# MINIMAL asBed, fromLine, toLine, chrom, chromStart, chromEnd #-}

-- | split a bed region into k consecutive subregions
splitBed :: BEDLike b => Int -> b -> [b]
splitBed k bed = map (uncurry (asBed chr)) . bins k $ (s, e)
  where
    chr = chrom bed
    s = chromStart bed
    e = chromEnd bed
{-# INLINE splitBed #-}

-- | split a bed region into consecutive fixed size subregions
splitBedBySize :: BEDLike b => Int -> b -> [b]
splitBedBySize k bed = map (uncurry (asBed chr)) . binBySize k $ (s, e)
  where
    chr = chrom bed
    s = chromStart bed
    e = chromEnd bed
{-# INLINE splitBedBySize #-}

-- | a type to imply that underlying data structure is sorted
newtype Sorted b = Sorted {fromSorted :: b}

compareBed :: (BEDLike b1, BEDLike b2) => b1 -> b2 -> Ordering
compareBed x y = compare x' y'
  where
    x' = (chrom x, chromStart x, chromEnd x)
    y' = (chrom y, chromStart y, chromEnd y)
{-# INLINE compareBed #-}

-- | sort BED, first by chromosome (alphabetical order), then by chromStart, last by chromEnd
sortBed :: BEDLike b => V.Vector b -> Sorted (V.Vector b)
sortBed beds = Sorted $ runST $ do
    v <- V.thaw beds
    I.sortBy compareBed v
    V.unsafeFreeze v
{-# INLINE sortBed #-}

mergeBed :: (BEDLike b, Monad m) => V.Vector b -> Source m b
mergeBed = mergeBed' . sortBed
{-# INLINE mergeBed #-}

mergeBed' :: (BEDLike b, Monad m) => Sorted (V.Vector b) -> Source m b
mergeBed' (Sorted beds) = source (V.head beds) 1
  where
    source !c !i
        | i >= n = yield c
        | otherwise = do
            let chr = chrom c
                s = chromStart c
                e = chromEnd c
                bed = beds V.! i
                chr' = chrom bed
                s' = chromStart bed
                e' = chromEnd bed
            case () of
                _ | chr /= chr' || s' > e+1 -> yield c >> source bed (i+1)
                  | e' <= e -> source c (i+1)
                  | otherwise -> source (asBed chr s e') (i+1)
    n = V.length beds
{-# INLINE mergeBed' #-}


-- * BED6 format

-- | BED6 format, as described in http://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
data BED = BED
    { _chrom :: !B.ByteString
    , _chromStart :: {-# UNPACK #-} !Int
    , _chromEnd :: {-# UNPACK #-} !Int
    , _name :: !(Maybe B.ByteString)
    , _score :: !(Maybe Double)
    , _strand :: !(Maybe Bool)  -- ^ True: "+", False: "-"
    } deriving (Eq, Show, Generic)

instance Default BED where
    def = BED
        { _chrom = ""
        , _chromStart = 0
        , _chromEnd = 0
        , _name = Nothing
        , _score = Nothing
        , _strand = Nothing
        }

instance BEDLike BED where
    asBed chr s e = BED chr s e Nothing Nothing Nothing

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

-- | retreive sequences
fetchSeq :: BioSeq DNA a => Genome -> Conduit BED IO (DNA a)
fetchSeq g = do gH <- lift $ gHOpen g
                table <- lift $ readIndex gH
                conduitWith gH table
                lift $ gHClose gH
  where
    conduitWith h index' = do 
        bed <- await
        case bed of
            Just (BED chr start end _ _ isForward) -> do 
                dna <- lift $ getSeq h index' (chr, start, end)
                case isForward of
                    Just False -> yield $ rc dna
                    _ -> yield dna
                conduitWith h index'
            _ -> return ()
{-# INLINE fetchSeq #-}


-- * BED3 format

data BED3 = BED3 !B.ByteString !Int !Int deriving (Eq, Show, Generic)

instance Default BED3 where
    def = BED3 "" 0 0

instance BEDLike BED3 where
    asBed = BED3

    fromLine l = case B.split '\t' l of
                    (a:b:c:_) -> BED3 a (readInt b) $ readInt c
                    _ -> error "Read BED fail: Incorrect number of fields"
    {-# INLINE fromLine #-}
    
    toLine (BED3 f1 f2 f3) = B.intercalate "\t" [f1, (B.pack.show) f2, (B.pack.show) f3]
    {-# INLINE toLine #-}

    chrom (BED3 f1 _ _) = f1
    chromStart (BED3 _ f2 _) = f2
    chromEnd (BED3 _ _ f3) = f3
