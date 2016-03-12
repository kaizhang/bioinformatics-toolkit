{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
module Bio.Motif.Alignment
    ( alignment
    , alignmentBy
    , linPenal
    , quadPenal
    , cubPenal
    , expPenal
    , l1
    , l2
    , l3
    , lInf
    , AlignFn
    , CombineFn
    ) where

import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.Unboxed as M
import Statistics.Sample (mean)

import Bio.Motif
import Bio.Utils.Functions

-- | penalty function takes the number of gaps and matched positions as input,
-- return penalty value
type PenalFn = Int -> Int -> Double

type DistanceFn = forall v. (G.Vector v Double, G.Vector v (Double, Double))
               => v Double -> v Double -> Double

type AlignFn = PWM
            -> PWM
            -> (Double, (Bool, Int))  -- ^ (distance, (on same direction,
                                      -- position w.r.t. the first pwm))

-- | combine distances from different positions of alignment
type CombineFn = U.Vector Double -> Double

alignment :: AlignFn
alignment = alignmentBy jsd (expPenal 0.05) l1

-- | linear penalty
linPenal :: Double -> PenalFn
linPenal x nGap nMatch = fromIntegral nGap * x / fromIntegral nMatch
{-# INLINE linPenal #-}

-- | quadratic penalty
quadPenal :: Double -> PenalFn
quadPenal x nGap nMatch = fromIntegral (nGap ^ (2 :: Int)) * x / fromIntegral nMatch
{-# INLINE quadPenal #-}

-- | cubic penalty
cubPenal :: Double -> PenalFn
cubPenal x nGap nMatch = fromIntegral (nGap ^ (3 :: Int)) * x / fromIntegral nMatch
{-# INLINE cubPenal #-}

-- | exponentail penalty
expPenal :: Double -> PenalFn
expPenal x nGap nMatch = fromIntegral (2^nGap :: Int) * x / fromIntegral nMatch
{-# INLINE expPenal #-}

l1 :: CombineFn
l1 = mean
{-# INLINE l1 #-}

l2 :: CombineFn
l2 = sqrt . mean . U.map (**2)
{-# INLINE l2 #-}

l3 :: CombineFn
l3 = (**(1/3)) . mean . U.map (**3)
{-# INLINE l3 #-}

lInf :: CombineFn
lInf = U.maximum
{-# INLINE lInf #-}

-- internal gaps are not allowed, larger score means larger distance, so the smaller the better
alignmentBy :: DistanceFn  -- ^ compute the distance between two aligned pwms
            -> PenalFn     -- ^ gap penalty
            -> CombineFn
            -> AlignFn
alignmentBy fn pFn combFn m1 m2
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
        | opti U.! i >= min' = (min',i')
        | d < min' = loop opti (d,i) a xs (i+1)
        | otherwise = loop opti (min',i') a xs (i+1)
      where
        d = combFn sc + pFn nGap nMatch
        sc = U.fromList $ zipWith fn a b
        nMatch = U.length sc
        nGap = n1 + n2 - 2 * nMatch
    loop _ acc _ _ _ = acc

    opti1 = optimalSc n1 n2
    opti2 = optimalSc n2 n1

    optimalSc x y = U.fromList $ scanr1 f $ go 0
      where
        f v min' = min v min'
        go i | nM == 0 = []
             | otherwise = pFn nG nM : go (i+1)
          where
            nM = min x $ y - i
            nG = i + abs (x - (y-i))

    s1 = M.toRows . _mat $ m1
    s2 = M.toRows . _mat $ m2
    s2' = M.toRows . _mat $ m2'
    m2' = rcPWM m2
    n1 = length s1
    n2 = length s2
{-# INLINE alignmentBy #-}
