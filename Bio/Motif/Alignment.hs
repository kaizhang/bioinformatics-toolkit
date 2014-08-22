{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE FlexibleContexts #-}

module Bio.Motif.Alignment where

import Bio.Motif
import Data.Packed.Matrix
import Data.Function (on)
import qualified Data.Vector.Generic as G
import Data.Clustering.Hierarchical
import NLP.Scores

type PenalFn = Int -> Double

type DistanceFn = G.Vector v Double => v Double -> v Double -> Double

-- | calculate distance between PWMs
distanceBy :: DistanceFn -> PWM -> PWM -> Double
distanceBy f (PWM _ m1) (PWM _ m2) = mean $ zipWith f x y
  where (x, y) = (toRows m1, toRows m2)
{-# INLINE distanceBy #-}

-- Calculate distance by Kullback-Leibler divergence
deltaKL :: PWM -> PWM -> Double
deltaKL = distanceBy (kullbackLeibler `on` G.toList)
{-# INLINE deltaKL #-}

-- Calculate distance by Jensen-Shannon divergence
deltaJS :: PWM -> PWM -> Double
deltaJS = distanceBy (jensenShannon `on` G.toList)
{-# INLINE deltaJS #-}

alignment :: PWM -> PWM -> (Double, (PWM, PWM, Int))
alignment = alignmentBy (jensenShannon `on` G.toList) quadPenal

linPenal :: PenalFn
linPenal n = fromIntegral n * 0.01

quadPenal :: PenalFn
quadPenal n = fromIntegral (n ^ (2 :: Int)) * 0.01

cubePenal :: PenalFn
cubePenal n = fromIntegral (n ^ (3 :: Int)) * 0.01

expPenal :: PenalFn
expPenal n = fromIntegral (2^n :: Int) * 0.01

-- internal gaps are not allowed, larger score means larger distance, so the smaller the better
alignmentBy :: DistanceFn -> PenalFn -> PWM -> PWM -> (Double, (PWM, PWM, Int))
alignmentBy fn pFn m1 m2 = ( fst result, (m1, m2, snd result) )
  where
    result = minimum $ zip (map (f s2 . flip drop s1) [0 .. n1-1]) [0 .. n1-1]
                    ++ zip (map (f s1 . flip drop s2) [1 .. n2-1]) [-1, -2 .. -n2+1]
    f a b = let xs = zipWith fn a b
                nGaps = n1 + n2 - 2 * length xs
            in mean xs + pFn nGaps
    s1 = toRows . _mat $ m1
    s2 = toRows . _mat $ m2
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
    -- shift motifs by 2
    δ (Motif _ x) (Motif _ y) = minimum $ map (deltaJS x) [shiftPWM y 0, shiftPWM y 1, shiftPWM y 2] ++ map (deltaJS y) [shiftPWM x 1, shiftPWM x 2]
    shiftPWM pwm i = subPWM i (n-i) pwm
      where n = rows $ _mat pwm
