{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}

module Bio.Motif.Alignment
    ( alignment
    , alignmentBy
    , mergePWM
    , buildTree
    , progressiveMerging
    ) where

import AI.Clustering.Hierarchical
import qualified Data.Vector.Generic as G
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Statistics.Matrix hiding (map)

import Bio.Motif
import Bio.Utils.Functions

-- | penalty function takes the gaps number as input, return penalty value
type PenalFn = Int -> Double

type DistanceFn = (G.Vector v Double, G.Vector v (Double, Double)) => v Double -> v Double -> Double

alignment :: PWM -> PWM -> (Double, (PWM, PWM, Int))
alignment = alignmentBy jsd quadPenal

-- | linear penalty
linPenal :: PenalFn
linPenal n = fromIntegral n * 0.3

-- | quadratic penalty
quadPenal :: PenalFn
quadPenal n = fromIntegral (n ^ (2 :: Int)) * 0.15

-- | cubic penalty
cubePenal :: PenalFn
cubePenal n = fromIntegral (n ^ (3 :: Int)) * 0.01

-- | exponentail penalty
expPenal :: PenalFn
expPenal n = fromIntegral (2^n :: Int) * 0.01

-- internal gaps are not allowed, larger score means larger distance, so the smaller the better
alignmentBy :: DistanceFn -> PenalFn -> PWM -> PWM -> (Double, (PWM, PWM, Int))
alignmentBy fn pFn m1 m2
    | fst forwardAlign <= fst reverseAlign = (fst forwardAlign, (m1, m2, snd forwardAlign))
    | otherwise = (fst reverseAlign, (m1, m2', snd reverseAlign))
  where
    forwardAlign | d1 < d2 = (d1,i1)
                 | otherwise = (d2,-i2)
      where
        (d1,i1) = loop (1/0,-1) s2 s1 0
        (d2,i2) = loop (1/0,-1) s1 s2 0
    reverseAlign | d1 < d2 = (d1,i1)
                 | otherwise = (d2,-i2)
      where
        (d1,i1) = loop (1/0,-1) s2' s1 0
        (d2,i2) = loop (1/0,-1) s1 s2' 0

    loop (min',i') a b@(_:xs) !i
        | currentBest >= min' = (min',i')
        | d < min' = loop (d,i) a xs (i+1)
        | otherwise = loop (min',i') a xs (i+1)
      where
        d = (G.sum sc + gapP) / fromIntegral (U.length sc + nGaps)
        currentBest = gapP / fromIntegral (n1 + n2)
        sc = U.fromList $ zipWith fn a b
        nGaps = n1 + n2 - 2 * U.length sc
        gapP = pFn nGaps
    loop (min',i') _ _ _ = (min',i')

    s1 = toRows . _mat $ m1
    s2 = toRows . _mat $ m2
    s2' = toRows . _mat $ m2'
    m2' = rcPWM m2
    n1 = length s1
    n2 = length s2
{-# INLINE alignmentBy #-}

mergePWM :: (PWM, PWM, Int) -> PWM
mergePWM (m1, m2, i) | i >= 0 = PWM Nothing (fromRows $ take i s1 ++ zipWith f (drop i s1) s2 ++ drop (n1 - i) s2)
                     | otherwise = PWM Nothing (fromRows $ take (-i) s2 ++ zipWith f (drop (-i) s2) s1 ++ drop (n2 + i) s1)
  where
    f = G.zipWith (\x y -> (x+y)/2)
    s1 = toRows . _mat $ m1
    s2 = toRows . _mat $ m2
    n1 = length s1
    n2 = length s2

progressiveMerging :: Dendrogram Motif -> PWM
progressiveMerging t = case t of
    Branch _ _ left right -> f (progressiveMerging left) $ progressiveMerging right
    Leaf a -> _pwm a
  where
    f a b = mergePWM $! snd $ alignment a b

-- | build a guide tree from a set of motifs
buildTree :: [Motif] -> Dendrogram Motif
buildTree motifs = hclust Average (V.fromList motifs) δ
  where
    δ (Motif _ x) (Motif _ y) = fst $ alignment x y
