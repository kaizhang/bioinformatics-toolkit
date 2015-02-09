{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE FlexibleContexts #-}

module Bio.Motif.Alignment where

import Bio.Motif
import Bio.Utils.Functions
import Data.Clustering.Hierarchical
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import Statistics.Matrix hiding (map)
import Statistics.Sample

-- | penalty function takes the gaps number as input, return penalty value
type PenalFn = Int -> Double

type DistanceFn = (G.Vector v Double, G.Vector v (Double, Double)) => v Double -> v Double -> Double

-- | calculate distance between PWMs
distanceBy :: DistanceFn -> PWM -> PWM -> Double
distanceBy f (PWM _ m1) (PWM _ m2) = mean . U.fromList $ zipWith f x y
  where (x, y) = (toRows m1, toRows m2)
{-# INLINE distanceBy #-}

-- Calculate distance by Kullback-Leibler divergence
deltaKL :: PWM -> PWM -> Double
deltaKL = distanceBy kld
{-# INLINE deltaKL #-}

-- Calculate distance by Jensen-Shannon divergence
deltaJS :: PWM -> PWM -> Double
deltaJS = distanceBy jsd
{-# INLINE deltaJS #-}

alignment :: PWM -> PWM -> (Double, (PWM, PWM, Int))
alignment = alignmentBy jsd 0.5

-- | linear penalty
linPenal :: PenalFn
linPenal n = fromIntegral n * 0.01

-- | quadratic penalty
quadPenal :: PenalFn
quadPenal n = fromIntegral (n ^ (2 :: Int)) * 0.01

-- | cubic penalty
cubePenal :: PenalFn
cubePenal n = fromIntegral (n ^ (3 :: Int)) * 0.01

-- | exponentail penalty
expPenal :: PenalFn
expPenal n = fromIntegral (2^n :: Int) * 0.01

-- internal gaps are not allowed, larger score means larger distance, so the smaller the better
alignmentBy :: DistanceFn -> Double -> PWM -> PWM -> (Double, (PWM, PWM, Int))
alignmentBy fn pFn m1 m2 | fst forwardAlign <= fst reverseAlign = (fst forwardAlign, (m1, m2, snd forwardAlign))
                         | otherwise = (fst reverseAlign, (m1, m2', snd reverseAlign))
  where
    forwardAlign = minimum $ zip (map (f s2 . flip drop s1) [0 .. n1-1]) [0 .. n1-1]
                          ++ zip (map (f s1 . flip drop s2) [1 .. n2-1]) [-1, -2 .. -n2+1]
    reverseAlign = minimum $ zip (map (f s2' . flip drop s1) [0 .. n1-1]) [0 .. n1-1]
                          ++ zip (map (f s1 . flip drop s2') [1 .. n2-1]) [-1, -2 .. -n2+1]
    f a b = let xs = U.fromList $ zipWith fn a b
                nGaps = n1 + n2 - 2 * U.length xs
            in (G.sum xs + pFn * fromIntegral nGaps) / fromIntegral (U.length xs + nGaps)
    s1 = toRows . _mat $ m1
    s2 = toRows . _mat $ m2
    s2' = toRows . _mat $ m2'
    m2' = rcPWM m2
    n1 = length s1
    n2 = length s2
{-# INLINE alignmentBy #-}

merge :: (PWM, PWM, Int) -> PWM
merge (m1, m2, i) | i >= 0 = PWM Nothing (fromRows $ take i s1 ++ zipWith f (drop i s1) s2 ++ drop (n1 - i) s2)
                  | otherwise = PWM Nothing (fromRows $ take (-i) s2 ++ zipWith f (drop (-i) s2) s1 ++ drop (n2 + i) s1)
  where
    f = G.zipWith (\x y -> (x+y)/2)
    s1 = toRows . _mat $ m1
    s2 = toRows . _mat $ m2
    n1 = length s1
    n2 = length s2

progressiveMerging :: Dendrogram Motif -> PWM
progressiveMerging t = case t of
    Branch _ left right -> f (progressiveMerging left) $ progressiveMerging right
    Leaf a -> _pwm a
  where
    f a b = merge $! snd $ alignment a b

-- | build a guide tree from a set of motifs
buildTree :: [Motif] -> Dendrogram Motif
buildTree motifs = dendrogram UPGMA motifs δ
  where
    δ (Motif _ x) (Motif _ y) = fst $ alignment x y
