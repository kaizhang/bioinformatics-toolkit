{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
module Bio.Motif.Alignment
    ( alignment
    , alignmentBy
    , linPenal
    , quadPenal
    , cubePenal
    , expPenal
    ) where

import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.Unboxed as M

import Bio.Motif
import Bio.Utils.Functions

-- | penalty function takes the gaps number as input, return penalty value
type PenalFn = Int -> Double

type DistanceFn = forall v. (G.Vector v Double, G.Vector v (Double, Double))
               => v Double -> v Double -> Double

alignment :: PWM -> PWM -> (Double, (Bool, Int))
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
alignmentBy :: DistanceFn  -- ^ compute the distance between two aligned pwms
            -> PenalFn     -- ^ gap penalty
            -> PWM
            -> PWM
            -> (Double, (Bool, Int))  -- ^ (distance, (on same direction,
                                      -- position w.r.t. the first pwm))
alignmentBy fn pFn m1 m2
    | fst forwardAlign <= fst reverseAlign =
        (fst forwardAlign, (True, snd forwardAlign))
    | otherwise = (fst reverseAlign, (False, snd reverseAlign))
  where
    forwardAlign | d1 < d2 = (d1,i1)
                 | otherwise = (d2,-i2)
      where
        (d1,i1) = loop opti2 (1/0,-1) s2 s1 0
        (d2,i2) = loop opti1 (1/0,-1) s1 s2 0
    reverseAlign | d1 < d2 = (d1,i1)
                 | otherwise = (d2,-i2)
      where
        (d1,i1) = loop opti2 (1/0,-1) s2' s1 0
        (d2,i2) = loop opti1 (1/0,-1) s1 s2' 0

    loop opti (min',i') a b@(_:xs) !i
        | currentBest >= min' = (min',i')
        | d < min' = loop opti (d,i) a xs (i+1)
        | otherwise = loop opti (min',i') a xs (i+1)
      where
        d = (G.sum sc + gapP) / fromIntegral (U.length sc + nGaps)
        currentBest = opti U.! i
        sc = U.fromList $ zipWith fn a b
        nGaps = n1 + n2 - 2 * U.length sc
        gapP = pFn nGaps
    loop _ (min',i') _ _ _ = (min',i')

    opti1 = optimalSc n1 n2
    opti2 = optimalSc n2 n1

    optimalSc x y = U.fromList $ scanr1 f $ go 0
      where
        f v min' = min v min'
        go i | a == 0 = []
             | otherwise = pFn b / fromIntegral (a + b) : go (i+1)
          where
            a = min x (y-i)
            b = i + abs (x - (y-i))

    s1 = M.toRows . _mat $ m1
    s2 = M.toRows . _mat $ m2
    s2' = M.toRows . _mat $ m2'
    m2' = rcPWM m2
    n1 = length s1
    n2 = length s2
{-# INLINE alignmentBy #-}
