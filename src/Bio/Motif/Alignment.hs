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
    , mergePWM
    , mergePWMWeighted
    , buildTree
    , mergeTree
    , mergeTreeWeighted
    ) where

import Control.Arrow (first)
import AI.Clustering.Hierarchical hiding (size)
import Data.List
import Data.Ord
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector.Generic as G
import qualified Data.Vector as V
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

mergePWM :: (PWM, PWM, Int) -> PWM
mergePWM (m1, m2, i) | i >= 0 = PWM Nothing (M.fromRows $ take i s1 ++ zipWith f (drop i s1) s2 ++ drop (n1 - i) s2)
                     | otherwise = PWM Nothing (M.fromRows $ take (-i) s2 ++ zipWith f (drop (-i) s2) s1 ++ drop (n2 + i) s1)
  where
    f = G.zipWith (\x y -> (x+y)/2)
    s1 = M.toRows . _mat $ m1
    s2 = M.toRows . _mat $ m2
    n1 = length s1
    n2 = length s2

mergePWMWeighted :: (PWM, [Int])  -- ^ pwm and weights at each position
                 -> (PWM, [Int])
                 -> Int                  -- ^ shift
                 -> (PWM, [Int])
mergePWMWeighted (pwm1, w1) (pwm2, w2) shift
    | shift >= 0 = toPWM $ take shift rows1 ++
                   zipWith f (drop shift rows1) rows2 ++ drop (n1 - shift) rows2
    | otherwise = toPWM $ take (-shift) rows2 ++
                  zipWith f (drop (-shift) rows2) rows1 ++ drop (n2 + shift) rows1
  where
    toPWM xs = let (rs, ws) = unzip xs
               in (PWM Nothing $ M.fromRows rs, ws)
    f (xs, a) (ys, b) = (G.zipWith (\x y -> (fromIntegral a * x + fromIntegral b * y) / fromIntegral (a + b)) xs ys, a + b)
    rows1 = zip (M.toRows $ _mat pwm1) w1
    rows2 = zip (M.toRows $ _mat pwm2) w2
    n1 = M.rows $ _mat pwm1
    n2 = M.rows $ _mat pwm2

mergeTree :: Dendrogram Motif -> PWM
mergeTree t = case t of
    Branch _ _ left right -> f (mergeTree left) $ mergeTree right
    Leaf a -> _pwm a
  where
    f a b | isSame = mergePWM (a, b, i)
          | otherwise = mergePWM (a, rcPWM b, i)
      where (_, (isSame, i)) = alignment a b

mergeTreeWeighted :: Dendrogram Motif -> (PWM, [Int])
mergeTreeWeighted t = case t of
    Branch _ _ left right -> f (mergeTreeWeighted left) $ mergeTreeWeighted right
    Leaf a -> (_pwm a, replicate (size $ _pwm a) 1)
  where
    f (a,w1) (b,w2)
        | isSame = mergePWMWeighted (a,w1) (b,w2) i
        | otherwise = mergePWMWeighted (a,w1) (rcPWM b, reverse w2) i
      where (_, (isSame, i)) = alignment a b
{-# INLINE mergeTreeWeighted #-}

iterativeMerge :: Double
               -> [Motif]   -- ^ Motifs to be merged. Motifs must have unique name.
               -> [(Motif, [Int])]
iterativeMerge th motifs = iter motifs'
  where
    motifs' = map (\x -> (x, replicate 1 $ size $ _pwm x)) motifs
    iter xs | d < th = iter $ first (Motif newName) merged :
                       filter (\x -> _name (fst x) /= _name x1 &&
                                     _name (fst x) /= _name x2 ) xs
            | otherwise = xs
      where
        ((d, (isSame, pos)), ((x1, w1), (x2, w2))) =
            minimumBy (comparing (fst . fst)) $
            map (\(a, b) -> (alignment (_pwm $ fst a) $ _pwm $ fst b, (a, b))) $ comb xs
        merged | isSame = mergePWMWeighted (_pwm x1, w1) (_pwm x2, w2) pos
               | otherwise = mergePWMWeighted (_pwm x1, w1) (rcPWM $ _pwm x2, reverse w2) pos
        newName = _name x1 `B.append` "+" `B.append` _name x2
    comb (y:ys) = zip (repeat y) ys ++ comb ys
    comb _ = []
{-# INLINE iterativeMerge #-}

-- | build a guide tree from a set of motifs
buildTree :: [Motif] -> Dendrogram Motif
buildTree motifs = hclust Average (V.fromList motifs) δ
  where
    δ (Motif _ x) (Motif _ y) = fst $ alignment x y

cutTreeBy :: Double   -- ^ start
          -> Double   -- ^ step
          -> ([Dendrogram a] -> Bool) -> Dendrogram a -> [Dendrogram a]
cutTreeBy start step fn tree = go start
  where
    go x | fn clusters = clusters
         | x - step > 0 = go $ x - step
         | otherwise = clusters
      where clusters = cutAt tree x
