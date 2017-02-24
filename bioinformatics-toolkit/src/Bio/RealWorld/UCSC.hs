{-# LANGUAGE OverloadedStrings #-}

module Bio.RealWorld.UCSC
    ( UCSCGene(..)
    , getTSS
    , getJunction
    , readUCSCGenes
    , readUCSCGenes'
    ) where

import qualified Data.ByteString.Char8 as B
import Conduit
import qualified Data.Vector.Unboxed as U
import System.IO

import Bio.RealWorld.ID
import Bio.Utils.Misc (readInt)

data UCSCGene = UCSCGene
    { _geneName :: !B.ByteString
    , _chrom :: !B.ByteString
    , _strand :: !Bool
    , _transcript :: !(Int, Int)
    , _cds :: !(Int, Int)
    , _exons :: !(U.Vector (Int, Int))
    , _introns :: !(U.Vector (Int, Int))
    , _proteinId :: !UniprotID
    , _alignId :: !UCSCID
    } deriving (Show)

-- | get Transcription Start Site
getTSS :: UCSCGene -> (B.ByteString, Int)
getTSS g = (_chrom g, fst $ _transcript g)

-- | get exon-intron junctions
getJunction :: UCSCGene -> (B.ByteString, U.Vector Int)
getJunction g = (_chrom g, U.map fst $ _introns g)

-- | read genes from UCSC "knownGenes.tsv"
readUCSCGenes :: FilePath -> Source IO UCSCGene
readUCSCGenes fl = do
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
{-# INLINE readUCSCGenes #-}

readUCSCGenes' :: FilePath -> IO [UCSCGene]
readUCSCGenes' fl = readUCSCGenes fl $$ sinkList
{-# INLINE readUCSCGenes' #-}

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
    in UCSCGene f1 f2 str trans cds exons introns (UniprotID f11) (UCSCID f12)
{-# INLINE readGeneFromLine #-}
