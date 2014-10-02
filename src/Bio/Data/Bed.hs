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

module Bio.Data.Bed (
      BED(..)
    , fetchSeq
    , readBED
    , readBED'
    , writeBED 
) where

import qualified Data.ByteString.Char8 as B
import Bio.Seq
import Bio.Seq.IO
import Bio.Utils.Misc (readInt, readDouble)
import Data.Conduit
import qualified Data.Conduit.List as CL
import Data.Maybe
import Control.Monad.State.Strict
import Data.Default.Generics
import GHC.Generics
import System.IO

-- | the type for BED format, as described in http://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
data BED = BED
    { _chrom :: !B.ByteString
    , _chromStart :: {-# UNPACK #-} !Int
    , _chromEnd :: {-# UNPACK #-} !Int
    , _name :: !(Maybe B.ByteString)
    , _score :: !(Maybe Double)
    , _strand :: !(Maybe Bool)  -- ^ True: "+", False: "-"
    } deriving (Read, Show, Generic)


instance Default BED

readBED :: FilePath -> Source IO BED
readBED fl = do handle <- liftIO $ openFile fl ReadMode
                loop handle
  where
    loop h = do eof <- liftIO $ hIsEOF h
                if eof 
                   then liftIO $ hClose h
                   else do
                       line <- liftIO $ B.hGetLine h
                       yield $ fromLine line
                       loop h
{-# INLINE readBED #-}

-- | non-streaming version
readBED' :: FilePath -> IO [BED]
readBED' fl = readBED fl $$ CL.consume

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

writeBED :: FilePath -> [BED] -> IO ()
writeBED fl beds = withFile fl WriteMode $ \h -> mapM_ (B.hPutStrLn h.toLine) beds
{-# INLINE writeBED #-}

fromLine :: B.ByteString -> BED
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

toLine :: BED -> B.ByteString
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
