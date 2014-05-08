{-# LANGUAGE OverloadedStrings, BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}

module Bio.ChIPSeq.ChIP
    ( naiveCC
    , naiveCC'
    ) where

import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import Data.List
import Bio.Utils.Bed
import Control.Parallel.Strategies
import Control.Lens ((^.))

fromBED :: [BED] -> [(B.ByteString, (S.HashSet Int, S.HashSet Int))]
{-# INLINE fromBED #-}
fromBED = map toSet.M.toList.M.fromListWith f.map parseLine
    where
        parseLine :: BED -> (B.ByteString, ([Int], [Int]))
        parseLine x = case x^.strand of
            Just True -> (x^.chrom, ([x^.chromStart], []))
            Just False -> (x^.chrom, ([], [x^.chromEnd]))
            _ -> error "Unknown Strand!"
        f (a,b) (a',b') = (a++a', b++b')
        toSet (chr, (forwd, rev)) = (chr, (S.fromList forwd, S.fromList rev))

-- | fast relative cross-correlation with smoothing
apprxCorr :: S.HashSet Int -> S.HashSet Int -> Int -> Int -> Int
{-# INLINE apprxCorr #-}
apprxCorr forwd rev smooth d = S.foldl' f 0 rev
    where
        f :: Int -> Int -> Int
        f !acc r | any (`S.member` forwd) [r-d-smooth..r-d+smooth] = acc + 1
                 | otherwise = acc

naiveCC :: [BED] -> [Int] -> Int
naiveCC = naiveCCWithSmooth 0

naiveCC' :: [BED] -> [Int] -> Int
naiveCC' = naiveCCWithSmooth 4

-- | calcuate cross corrlation with different shifts
naiveCCWithSmooth :: Int -> [BED] -> [Int] -> Int
{-# INLINE naiveCCWithSmooth #-}
naiveCCWithSmooth smooth input range = maxCC.fromBED $ input
    where
        cc :: [(B.ByteString, (S.HashSet Int, S.HashSet Int))] -> Int -> Int
        cc xs d = foldl' (+) 0 . map ((\(forwd, rev) -> apprxCorr forwd rev smooth d).snd) $ xs
        maxCC xs = snd.maximum $ zip (parMap rpar (cc xs) range) range

