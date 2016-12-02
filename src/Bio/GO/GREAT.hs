{-# LANGUAGE BangPatterns #-}
module Bio.GO.GREAT
    ( AssocRule(..)
    , defaultAssocRule
    , getRegulatoryDomains
    , read3DContact
    , get3DRegulatoryDomains
    ) where

import           Conduit
import           Control.Monad                (forM_)
import qualified Data.ByteString.Char8        as B
import           Data.Function                (on)
import qualified Data.IntervalMap             as IM
import           Data.List                    (foldl', sortBy)

import           Bio.Data.Bed
import           Bio.Utils.Misc               (readInt)

-- | How to associate genomic regions with genes
data AssocRule =
    BasalPlusExtension Int Int Int  -- ^ Upstream, downstream and extension
  | TwoNearest                      -- ^ Not implemented
  | OneNearest                      -- ^ Not implemented

-- | The default rule is defined as -5k ~ +1k around TSSs,
-- plus up to 1000k extension.
defaultAssocRule :: AssocRule
defaultAssocRule = BasalPlusExtension 5000 1000 1000000

-- | A Gene consists of the chromosome name, TSS, strandness and an associated value.
type Gene a = ((B.ByteString, Int, Bool), a)

-- | Given a gene list and the rule, compute the rulatory domain for each gene
getRegulatoryDomains :: AssocRule -> [Gene a] -> [(BED3, a)]
getRegulatoryDomains (BasalPlusExtension up dw ext) genes = (extendTail r ext, a) : rs
  where
    (rs, Just (r,a)) = foldl' f ([], Nothing) $ sortBy (compareBed `on` fst) basal
    f (acc, Nothing) (b,x) = (acc, Just (extendHead b ext, x))
    f (acc, Just (b', x')) (b,x)
        | chrom b' /= chrom b = ( (extendTail b' ext, x') : acc
                                , Just (extendHead b ext, x) )
        | chromEnd b' >= chromStart b = ((b',x') : acc, Just (b,x))
        | otherwise = let ext' = min ext $ (chromStart b - chromEnd b') `div` 2
                      in ((extendTail b' ext', x') : acc, Just (extendHead b ext', x))
    extendHead (BED3 chr s e) l | s - l >= 0 = BED3 chr (s-l) e
                                | otherwise = BED3 chr 0 e
    extendTail (BED3 chr s e) l = BED3 chr s (e+l)
    basal = flip map genes $ \((chr, tss, str), x) ->
        if str then (BED3 chr (tss - up) (tss + dw), x)
               else (BED3 chr (tss - dw) (tss + up), x)
getRegulatoryDomains _ _ = undefined
{-# INLINE getRegulatoryDomains #-}

-- | Read 3D contacts from a file, where each line contains 6 fields separated
-- by Tabs corresponding to 2 interacting loci. Example:
-- chr1 [TAB] 11 [TAB] 100 [TAB] chr2 [TAB] 23 [TAB] 200
read3DContact :: FilePath -> Source (ResourceT IO) (BED3, BED3)
read3DContact input = sourceFileBS input =$= linesUnboundedAsciiC =$= mapC f
  where
    f x = let (chr1:s1:e1:chr2:s2:e2:_) = B.split '\t' x
          in (BED3 chr1 (readInt s1) (readInt e1), BED3 chr2 (readInt s2) (readInt e2))

-- | 3D regulatory domains of a gene include the gene's promoter and regions
-- that contact with the gene's promoter.
get3DRegulatoryDomains :: Monad m
                       => [Gene a]   -- ^ Genes
                       -> Int        -- ^ Upstream
                       -> Int        -- ^ Downstream
                       -> Conduit (BED3, BED3) m (BED3, a)
get3DRegulatoryDomains genes up dw = concatRegions >> promoters
  where
    promoters = forM_ genes $ \((chr, tss, str), x) ->
        if str then yield (BED3 chr (gateZero $ tss - up) (tss + dw), x)
               else yield (BED3 chr (gateZero $ tss - dw) (tss + up), x)
    concatRegions = concatMapC $ \(locA, locB) ->
        map ((,) locB) (IM.elems $ intersecting basal locA) ++
        map ((,) locA) (IM.elems $ intersecting basal locB)
    basal = bedToTree undefined $ flip map genes $ \((chr, tss, str), x) ->
        if str then (BED3 chr (tss - up) (tss + dw), x)
               else (BED3 chr (tss - dw) (tss + up), x)
    gateZero x | x < 0 = 0
               | otherwise = x
{-# INLINE get3DRegulatoryDomains #-}
