{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
module Bio.Motif.Merge
    ( mergePWM
    , mergePWMWeighted
    , dilute
    , trim
    , mergeTree
    , mergeTreeWeighted
    , iterativeMerge
    , buildTree
    )where

import           AI.Clustering.Hierarchical    hiding (size)
import Control.Arrow (first)
import           Control.Monad                 (forM_, when)
import           Control.Monad.ST              (runST, ST)
import qualified Data.ByteString.Char8         as B
import           Data.List                     (dropWhileEnd)
import qualified Data.Matrix.Symmetric.Mutable as MSU
import qualified Data.Matrix.Unboxed           as MU
import           Data.Maybe
import qualified Data.Vector                   as V
import qualified Data.Vector.Mutable           as VM
import qualified Data.Vector.Unboxed           as U

import           Bio.Motif
import           Bio.Motif.Alignment
import           Bio.Utils.Functions           (kld)

mergePWM :: (PWM, PWM, Int) -> PWM
mergePWM (m1, m2, shift) | shift >= 0 = merge shift (_mat m1) $ _mat m2
                         | otherwise = merge (-shift) (_mat m2) $ _mat m1
  where
    merge s a b = PWM Nothing $ MU.fromRows $ loop 0
      where
        n1 = MU.rows a
        n2 = MU.rows b
        loop i | i' < 0 || (i < n1 && i' >= n2) = MU.takeRow a i : loop (i+1)
               | i < n1 && i' < n2 = f (MU.takeRow a i) (MU.takeRow b i') : loop (i+1)
               | i >= n1 && i' < n2 = MU.takeRow b i' : loop (i+1)
               | otherwise = []
          where
            i' = i - s
            f = U.zipWith (\x y -> (x+y)/2)

mergePWMWeighted :: (PWM, [Int])  -- ^ pwm and weights at each position
                 -> (PWM, [Int])
                 -> Int                  -- ^ shift
                 -> (PWM, [Int])
mergePWMWeighted m1 m2 shift
    | shift >= 0 = merge shift (first _mat m1) $ first _mat m2
    | otherwise = merge (-shift) (first _mat m2) $ first _mat m1

  where
    merge s (p1,w1) (p2,w2) = first (PWM Nothing . MU.fromRows) $ unzip $ loop 0
      where
        a = V.fromList $ zip (MU.toRows p1) w1
        b = V.fromList $ zip (MU.toRows p2) w2
        n1 = V.length a
        n2 = V.length b
        loop i | i' < 0 || (i < n1 && i' >= n2) = a V.! i : loop (i+1)
               | i < n1 && i' < n2 = f (a V.! i) (b V.! i') : loop (i+1)
               | i >= n1 && i' < n2 = b V.! i' : loop (i+1)
               | otherwise = []
          where
            i' = i - s
            f (xs, wx) (ys, wy) = (U.zipWith (\x y ->
                (fromIntegral wx * x + fromIntegral wy * y) /
                fromIntegral (wx + wy)) xs ys, wx + wy)

-- | dilute positions in a PWM that are associated with low weights
dilute :: (PWM, [Int]) -> PWM
dilute (pwm, ws) = PWM Nothing $ MU.fromRows $ zipWith f ws $ MU.toRows $ _mat pwm
  where
    f w r | w < n = let d = fromIntegral $ n - w
                    in U.map (\x -> (fromIntegral w * x + 0.25 * d) / fromIntegral n) r
          | otherwise = r
    n = maximum ws
{-# INLINE dilute #-}

trim :: Bkgd -> Double -> PWM -> PWM
trim (BG (a,c,g,t)) cutoff pwm = PWM Nothing $ MU.fromRows $ dropWhileEnd f $
    dropWhile f rs
  where
    f x = kld x bg < cutoff
    rs = MU.toRows $ _mat pwm
    bg = U.fromList [a,c,g,t]
{-# INLINE trim #-}

mergeTree :: AlignFn -> Dendrogram Motif -> PWM
mergeTree align t = case t of
    Branch _ _ left right -> f (mergeTree align left) $ mergeTree align right
    Leaf a -> _pwm a
  where
    f a b | isSame = mergePWM (a, b, i)
          | otherwise = mergePWM (a, rcPWM b, i)
      where (_, (isSame, i)) = align a b

mergeTreeWeighted :: AlignFn -> Dendrogram Motif -> (PWM, [Int])
mergeTreeWeighted align t = case t of
    Branch _ _ left right -> f (mergeTreeWeighted align left) $
        mergeTreeWeighted align right
    Leaf a -> (_pwm a, replicate (size $ _pwm a) 1)
  where
    f (a,w1) (b,w2)
        | isSame = mergePWMWeighted (a,w1) (b,w2) i
        | otherwise = mergePWMWeighted (a,w1) (rcPWM b, reverse w2) i
      where (_, (isSame, i)) = align a b
{-# INLINE mergeTreeWeighted #-}

iterativeMerge :: AlignFn
               -> Double    -- cutoff
               -> [Motif]   -- ^ Motifs to be merged. Motifs must have unique name.
               -> [([B.ByteString], PWM, [Int])]
iterativeMerge align th motifs = runST $ do
    motifs' <- V.unsafeThaw $ V.fromList $ flip map motifs $ \x ->
        Just ([_name x], _pwm x, replicate (size $ _pwm x) 1)

    let n = VM.length motifs'
        iter mat = do
            -- retrieve the minimum value
            ((i, j), (d, (isSame, pos))) <- loop ((-1,-1), (1/0, undefined)) 0 1
            if d < th
                then do
                    Just (nm1, pwm1, w1) <- VM.unsafeRead motifs' i
                    Just (nm2, pwm2, w2) <- VM.unsafeRead motifs' j
                    let merged = (nm1 ++ nm2, pwm', w')
                        (pwm',w') | isSame = mergePWMWeighted (pwm1, w1) (pwm2, w2) pos
                                  | otherwise = mergePWMWeighted (pwm1, w1)
                                              (rcPWM $ pwm2, reverse w2) pos

                    -- update
                    forM_ [0..n-1] $ \i' -> MSU.unsafeWrite mat (i',j) Nothing
                    VM.unsafeWrite motifs' i $ Just merged
                    VM.unsafeWrite motifs' j Nothing
                    forM_ [0..n-1] $ \j' -> when (i /= j') $ do
                        x <- VM.unsafeRead motifs' j'
                        case x of
                            Just (_, pwm2',_) -> do
                                let ali | i < j' = Just $ align pwm' pwm2'
                                        | otherwise = Just $ align pwm2' pwm'
                                MSU.unsafeWrite mat (i,j') ali
                            _ -> return ()
                    iter mat
                else return ()
          where
            loop ((i_min, j_min), d_min) i j
                | i >= n = return ((i_min, j_min), d_min)
                | j >= n = loop ((i_min, j_min), d_min) (i+1) (i+2)
                | otherwise = do
                    x <- MSU.unsafeRead mat (i,j)
                    case x of
                        Just d -> if fst d < fst d_min
                            then loop ((i,j), d) i (j+1)
                            else loop ((i_min, j_min), d_min) i (j+1)
                        _ -> loop ((i_min, j_min), d_min) i (j+1)

    -- initialization
    mat <- MSU.replicate (n,n) Nothing :: ST s (MSU.SymMMatrix VM.MVector s (Maybe (Double, (Bool, Int))))
    forM_ [0..n-1] $ \i -> forM_ [i+1 .. n-1] $ \j -> do
        Just (_, pwm1, _) <- VM.unsafeRead motifs' i
        Just (_, pwm2, _) <- VM.unsafeRead motifs' j
        MSU.unsafeWrite mat (i,j) $ Just $ align pwm1 pwm2

    iter mat
    results <- V.unsafeFreeze motifs'
    return $ V.toList $ V.map fromJust $ V.filter isJust results
{-# INLINE iterativeMerge #-}

-- | build a guide tree from a set of motifs
buildTree :: AlignFn -> [Motif] -> Dendrogram Motif
buildTree align motifs = hclust Average (V.fromList motifs) δ
  where
    δ (Motif _ x) (Motif _ y) = fst $ align x y

cutTreeBy :: Double   -- ^ start
          -> Double   -- ^ step
          -> ([Dendrogram a] -> Bool) -> Dendrogram a -> [Dendrogram a]
cutTreeBy start step fn tree = go start
  where
    go x | fn clusters = clusters
         | x - step > 0 = go $ x - step
         | otherwise = clusters
      where clusters = cutAt tree x
