{-# LANGUAGE OverloadedStrings #-}
--------------------------------------------------------------------------------
-- |
-- Module      :  $Header$
-- Description :  Search and download data from ENCODE project
-- Copyright   :  (c) Kai Zhang
-- License     :  MIT

-- Maintainer  :  kai@kzhang.org
-- Stability   :  experimental
-- Portability :  portable

-- resources from UCSC
--------------------------------------------------------------------------------

module Bio.RealWorld.UCSC
    ( UCSCGene(..)
    , readGenes
    , readGenes'
    ) where

import Control.Monad.IO.Class (liftIO)
import qualified Data.ByteString.Char8 as B
import Data.Conduit
import qualified Data.Conduit.List as CL
import qualified Data.Vector.Unboxed as U
import System.IO

import Bio.Utils.Misc (readInt)

data UCSCGene = UCSCGene
    { _geneName :: !B.ByteString
    , _chrom :: !B.ByteString
    , _strand :: !Bool
    , _transcript :: !(Int, Int)
    , _cds :: !(Int, Int)
    , _exons :: !(U.Vector (Int, Int))
    , _introns :: !(U.Vector (Int, Int))
    , _proteinId :: !B.ByteString
    , _alignId :: !B.ByteString
    } deriving (Show)

-- | read genes from UCSC "knownGenes.tsv"
readGenes :: FilePath -> Source IO UCSCGene
readGenes fl = do
    handle <- liftIO $ openFile fl ReadMode
    _ <- liftIO $ B.hGetLine handle   -- header
    loop handle
  where
    loop h = do
        eof <- liftIO $ hIsEOF h
        if eof
           then liftIO $ hClose h
           else do
               l <- liftIO $ B.hGetLine h
               yield $ readGeneFromLine l
               loop h
{-# INLINE readGenes #-}

readGenes' :: FilePath -> IO [UCSCGene]
readGenes' fl = readGenes fl $$ CL.consume
{-# INLINE readGenes' #-}

readGeneFromLine :: B.ByteString -> UCSCGene
readGeneFromLine xs =
    let [f1,f2,f3,f4,f5,f6,f7,_,f9,f10,f11,f12] = B.split '\t' xs
        str | f3 == "+" = True
            | otherwise = False
        trans = (readInt f4, readInt f5)
        cds = (readInt f6, readInt f7)
        exonStarts = map readInt . init . B.split ',' $ f9
        exonEnds = map readInt . init . B.split ',' $ f10
        exons = U.fromList $ zip exonStarts exonEnds
        introns = U.fromList $ zip exonEnds $ tail exonStarts
    in UCSCGene f1 f2 str trans cds exons introns f11 f12
{-# INLINE readGeneFromLine #-}

