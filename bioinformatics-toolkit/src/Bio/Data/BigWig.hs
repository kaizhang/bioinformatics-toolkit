module Bio.Data.BigWig
    ( module Data.BBI.BigWig
    , accumScores
    , normalizedScores
    ) where

import Data.BBI.BigWig
import Data.Conduit
import qualified Data.Conduit.List as CL
import Bio.Data.Bed

-- | for each BED region, return total scores of all wig records overlapped with it
accumScores :: BEDLike b => BWFile -> Conduit b IO Double
accumScores bwF = CL.mapM (helper bwF)
-}

-- | for each BED region, return normalized scores (divided by the length) of
-- all wig records overlapped with it
normalizedScores :: BEDLike b => BWFile -> Conduit b IO Double
normalizedScores bwF = CL.mapM $ \bed -> do
    x <- helper bwF bed
    return $ x / (fromIntegral . size) bed

helper :: BEDLike b => BWFile -> b -> IO Double
helper bwF bed = queryBWFile bwF (chr, start, end) $$ CL.fold f 0.0
  where
    f acc (_, s, e, v) = acc + fromIntegral (e - s) * v
    chr = chrom bed
    start = chromStart bed
    end = chromEnd bed
