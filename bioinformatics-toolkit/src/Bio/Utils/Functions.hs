{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}

module Bio.Utils.Functions (
      ihs
    , ihs'
    , scale
    , filterFDR
    , slideAverage
    , hyperquick
    , gaussianKDE
    , kld
    , jsd
    , binarySearch
    , binarySearchBy
    , binarySearchByBounds
    , quantileNormalization
    , quantileNormalization'
) where

import Data.Bits (shiftR)
import Data.List (foldl', groupBy)
import Data.Function (on)
import Data.Ord (comparing)
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix as M
import Statistics.Sample (meanVarianceUnb, mean)
import Statistics.Function (sortBy)
import Control.Parallel.Strategies (parMap, rseq)
import Statistics.Sample.KernelDensity (kde)

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

-- | given the p-values, filter data by controling FDR
filterFDR :: G.Vector v (a, Double)
          => Double  -- ^ desired FDR value
          -> v (a, Double)  -- ^ input data and pvalues
          -> v (a, Double)
filterFDR α xs = go n . sortBy (comparing snd) $ xs
  where
    go rank v | rank <= 0 = G.empty
              | snd (v `G.unsafeIndex` (rank-1)) <= fromIntegral rank * α / fromIntegral n = G.slice 0 rank v
              | otherwise = go (rank-1) v
    n = G.length xs
{-# INLINE filterFDR #-}

-- | Compute the sliding average for each entry in a vector
slideAverage :: (Fractional a, G.Vector v a)
             => Int              -- ^ size of HALF sliding window, 2 means a total
                                 -- window size is 5
             -> v a
             -> v a
slideAverage k xs = G.generate n $ \i -> go xs (max 0 (i-k)) (min (n-1) (i+k))
  where
    go v i j = let l = j - i + 1
               in G.foldl' (+) 0 (G.unsafeSlice i l v) / fromIntegral l
    n = G.length xs
{-# INLINE slideAverage #-}

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
    invJm _n _x _N _m = ( 1 - fromIntegral _x / fromIntegral (_m+1) ) /
                          ( 1 - fromIntegral (_n-1-_x) / fromIntegral (_N-1-_m) )
    e = 1e-20

-- | Assign weights to the points according to density estimation.
gaussianKDE :: Int  -- ^ number of mesh points used in KDE
            -> U.Vector Double -> (Double -> Double)
gaussianKDE n xs = \x -> 
    let i = binarySearch points x
        lo = points U.! (i - 1)
        lo_d = den U.! (i - 1) 
        hi = points U.! i
        hi_d = den U.! i
        hi_w = (x - lo) / (hi - lo)
        lo_w = 1 - hi_w
     in lo_w * lo_d + hi_w * hi_d
  where
    (points, den) = kde n xs

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

-- | O(log n). return the position of the first element that is >= query
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
                        EQ -> loop l (k-1)
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

-- | Columns are samples, rows are features / genes.
quantileNormalization :: M.Matrix Double -> M.Matrix Double
quantileNormalization mat = M.fromColumns $ map
    (fst . G.unzip . sortBy (comparing snd) . G.fromList . concatMap f .
    groupBy ((==) `on` (snd . snd)) . zip averages . G.toList) $
    M.toColumns srtMat
  where
    f [(a,(b,_))] = [(a,b)]
    f xs = let m = mean $ U.fromList $ fst $ unzip xs
           in map (\(_,(i,_)) -> (m, i)) xs
    srtMat :: M.Matrix (Int, Double)
    srtMat = M.fromColumns $ map (sortBy (comparing snd) . G.zip (G.enumFromN 0 n)) $
        M.toColumns mat
    averages = map (mean . snd . G.unzip) $ M.toRows srtMat
    n = M.rows mat
{-# INLINE quantileNormalization #-}

-- | Columns are samples, rows are features / genes.
quantileNormalization' :: M.Matrix Double -> M.Matrix Double
quantileNormalization' mat = M.fromColumns $ parMap rseq
    (fst . G.unzip . sortBy (comparing snd) . G.fromList . concatMap f .
    groupBy ((==) `on` (snd . snd)) . zip averages . G.toList) $
    M.toColumns srtMat
  where
    f [(a,(b,_))] = [(a,b)]
    f xs = let m = mean $ U.fromList $ fst $ unzip xs
           in map (\(_,(i,_)) -> (m, i)) xs
    srtMat :: M.Matrix (Int, Double)
    srtMat = M.fromColumns $ parMap rseq (sortBy (comparing snd) . G.zip (G.enumFromN 0 n)) $
        M.toColumns mat
    averages = parMap rseq (mean . snd . G.unzip) $ M.toRows srtMat
    n = M.rows mat
{-# INLINE quantileNormalization' #-}