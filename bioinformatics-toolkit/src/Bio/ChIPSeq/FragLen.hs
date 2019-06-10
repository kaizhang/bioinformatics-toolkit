{-# LANGUAGE BangPatterns #-}

module Bio.ChIPSeq.FragLen
    ( fragLength
    , naiveCCWithSmooth
    ) where

import           Bio.Data.Bed
import           Lens.Micro                ((^.))
import           Control.Parallel.Strategies (parMap, rpar)
import qualified Data.ByteString.Char8       as B
import qualified Data.HashMap.Strict         as M
import qualified Data.HashSet                as S
import           Data.List                   (foldl', maximumBy)
import           Data.Ord                    (comparing)

-- | estimate fragment length for a ChIP-seq experiment
fragLength :: (Int, Int) -> [BED] -> Int
fragLength (start, end) beds = fst $ maximumBy (comparing snd) $
    naiveCCWithSmooth 4 beds [start, start+2 .. end]
{-# INLINE fragLength #-}

-- sizeDistribution :: (Int, Int) -> [BED] -> (Int, Int)

fromBED :: [BED] -> [(B.ByteString, (S.HashSet Int, S.HashSet Int))]
fromBED = map toSet . M.toList . M.fromListWith f . map parseLine
    where
        parseLine :: BED -> (B.ByteString, ([Int], [Int]))
        parseLine x = case x^.strand of
            Just True  -> (x^.chrom, ([x^.chromStart], []))
            Just False -> (x^.chrom, ([], [x^.chromEnd]))
            _          -> error "Unknown Strand!"
        f (a,b) (a',b') = (a++a', b++b')
        toSet (chr, (forwd, rev)) = (chr, (S.fromList forwd, S.fromList rev))
{-# INLINE fromBED #-}

-- | fast relative cross-correlation with smoothing
apprxCorr :: S.HashSet Int -> S.HashSet Int -> Int -> Int -> Int
apprxCorr forwd rev smooth d = S.foldl' f 0 rev
    where
        f :: Int -> Int -> Int
        f !acc r | any (`S.member` forwd) [r-d-smooth..r-d+smooth] = acc + 1
                 | otherwise = acc
{-# INLINE apprxCorr #-}

-- | calcuate cross corrlation with different shifts
naiveCCWithSmooth :: Int -> [BED] -> [Int] -> [(Int, Int)]
naiveCCWithSmooth smooth input range = f . fromBED $ input
    where
        cc :: [(B.ByteString, (S.HashSet Int, S.HashSet Int))] -> Int -> Int
        cc xs d = foldl' (+) 0 . map ((\(forwd, rev) -> apprxCorr forwd rev smooth d).snd) $ xs
        f xs = zip range $ parMap rpar (cc xs) range
{-# INLINE naiveCCWithSmooth #-}
