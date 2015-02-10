{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
--------------------------------------------------------------------------------
-- |
-- Module      :  $Header$
-- Copyright   :  (c) 2014 Kai Zhang
-- License     :  MIT

-- Maintainer  :  kai@kzhang.org
-- Stability   :  experimental
-- Portability :  portable

-- some useful functions
--------------------------------------------------------------------------------

module Bio.Utils.Functions (
      ihs
    , ihs'
    , scale
    , hyperquick
    , kld
    , jsd
    , binarySearch
    , binarySearchBy
    , binarySearchByBounds
) where

import Data.Bits (shiftR)
import Data.List (foldl')
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import Statistics.Sample (meanVarianceUnb)

-- | inverse hyperbolic sine transformation
ihs :: Double  -- ^ θ, determine the shape of the function
    -> Double  -- ^ input
    -> Double
ihs !θ !x | θ == 0 = x
          | otherwise = log (θ * x + sqrt (θ * θ * x * x + 1)) / θ
{-# INLINE ihs #-}

-- | inverse hyperbolic sine transformation with θ = 1
ihs' :: Double -> Double
ihs' = ihs 1

-- | scale data to zero mean and 1 variance
scale :: G.Vector v Double => v Double -> v Double
scale xs = G.map (\x -> (x - m) / sqrt s) xs
  where
    (m,s) = meanVarianceUnb xs
{-# INLINE scale #-}

hyperquick :: Int -> Int -> Int -> Int -> Double
hyperquick x m _n _N = loop (m-2) s s (2*e)
  where
    loop !k !ak !bk !epsk
        | k < _N - (_n-x) - 1 && epsk > e =
            let ck = ak / bk
                k' = k + 1
                jjm = invJm _n x _N k'
                bk' = bk * jjm + 1
                ak' = ak * jjm
                espk' = fromIntegral (_N - (_n - x) - 1 - k') * (ck - ak' / bk')
            in loop k' ak' bk' espk'
        | otherwise = 1 - (ak / bk - epsk / 2)
    s = foldl' (\s' k -> 1 + s' * invJm _n x _N k) 1.0 [x..m-2]
    invJm _n x _N m = ( 1 - fromIntegral x / fromIntegral (m+1) ) /
                          ( 1 - fromIntegral (_n-1-x) / fromIntegral (_N-1-m) )
    e = 1e-20

{-
hyperquick' :: Int -> Int -> Int -> Int -> Double
hyperquick' x m _n _N = loop (m-2) s s (2*e)
  where
    loop !k !ak !bk !epsk
        | k < _N - (_n-x) - 1 && epsk > e =
            let ck = ak / bk
                k' = k + 1
                jjm = invJm _n x _N k'
                bk' = bk * jjm + 1
                ak' = ak * jjm
                espk' = fromIntegral (_N - (_n - x) - 1 - k') * (ck - ak' / bk')
            in loop k' ak' bk' espk'
        | otherwise = 1 - (ak / bk - epsk / 2)
    s = foldl' (\s' k -> 1 + s' * invJm _n x _N k) 1.0 [x..m-2]
    invJm _n x _N m = ( 1 - fromIntegral x / fromIntegral (m+1) ) /
                          ( 1 - fromIntegral (_n-1-x) / fromIntegral (_N-1-m) )
    e = 1e-200
    -}

-- | compute the Kullback-Leibler divergence between two valid (not check) probability distributions.
-- kl(X,Y) = \sum_i P(x_i) log_2(P(x_i)\/P(y_i)). 
kld :: (G.Vector v Double, G.Vector v (Double, Double)) => v Double -> v Double -> Double
kld xs ys | G.length xs /= G.length ys = error "Incompitable dimensions"
          | otherwise = G.foldl' f 0 . G.zip xs $ ys
  where
    f acc (x, y) | x == 0 = acc
                 | otherwise = acc + x * (logBase 2 x - logBase 2 y)
{-# SPECIALIZE kld :: U.Vector Double -> U.Vector Double -> Double #-}
{-# SPECIALIZE kld :: V.Vector Double -> V.Vector Double -> Double #-}

-- | Jensen-Shannon divergence: JS(X,Y) = 1\/2 KL(X,(X+Y)\/2) + 1\/2 KL(Y,(X+Y)\/2).
jsd :: (G.Vector v Double, G.Vector v (Double, Double)) => v Double -> v Double -> Double
jsd xs ys = 0.5 * kld xs zs + 0.5 * kld ys zs
  where zs = G.zipWith (\x y -> (x + y) / 2) xs ys
{-# SPECIALIZE jsd :: U.Vector Double -> U.Vector Double -> Double #-}
{-# SPECIALIZE jsd :: V.Vector Double -> V.Vector Double -> Double #-}

-- | O(log n). return the position of the first element that is greater than query
binarySearch :: (G.Vector v e, Ord e)
             => v e -> e -> Int
binarySearch vec e = binarySearchByBounds compare vec e 0 $ G.length vec - 1
{-# INLINE binarySearch #-}

binarySearchBy :: G.Vector v e
               => (e -> a -> Ordering) -> v e -> a -> Int
binarySearchBy cmp vec e = binarySearchByBounds cmp vec e 0 $ G.length vec - 1
{-# INLINE binarySearchBy #-}

binarySearchByBounds :: G.Vector v e
                     => (e -> a -> Ordering) -> v e -> a -> Int -> Int -> Int
binarySearchByBounds cmp vec e = loop
  where
    loop !l !u
        | u < l = l
        | otherwise = case cmp (vec G.! k) e of
                        LT -> loop (k+1) u
                        EQ -> k
                        GT -> loop l (k-1)
      where k = u + l `shiftR` 1
{-# INLINE binarySearchByBounds #-}

{-
empiricalCDF :: G.Vector v Double => v Double -> v Double
empiricalCDF xs = runST $ do
    let n = G.length xs
        indices = groupBy ( (==) `on` ((xs G.!).snd) ) $ zip [1.0..] $ sortBy (compare `on` (xs G.!)) [0..n-1]
        updates mv (v,ys) = mapM_ (flip (GM.unsafeWrite mv) v.snd) ys
    xs' <- G.thaw xs
    mapM_ (updates xs'. ((flip (/) (fromIntegral n).fst.last) &&& id)) indices
    G.unsafeFreeze xs'
{-# INLINE empiricalCDF #-}

cdf :: G.Vector v Double => v Double -> v Double
cdf xs = let f = empiricalCDF xs
             n = fromIntegral $ G.length xs
             δ = 1 / (4 * (n**0.25) * sqrt (pi * log n))
         in G.map (\ x -> case () of
             _ | x < δ -> δ
               | x > 1 - δ -> 1 - δ
               | otherwise -> x) f
{-# INLINE cdf #-}
-}
