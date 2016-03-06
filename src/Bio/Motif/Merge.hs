{-# LANGUAGE OverloadedStrings #-}
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
import           Control.Monad                 (forM_, when)
import           Control.Monad.ST              (runST)
import qualified Data.ByteString.Char8         as B
import           Data.List                     (dropWhileEnd)
import qualified Data.Matrix.Symmetric         as MS
import qualified Data.Matrix.Symmetric.Mutable as MSU
import qualified Data.Matrix.Unboxed           as MU
import           Data.Maybe
import           Data.Ord                      (comparing)
import qualified Data.Vector                   as V
import qualified Data.Vector.Mutable           as VM
import qualified Data.Vector.Unboxed           as U

import           Bio.Motif
import           Bio.Motif.Alignment
import           Bio.Utils.Functions           (kld)

mergePWM :: (PWM, PWM, Int) -> PWM
mergePWM (m1, m2, i)
    | i >= 0 = PWM Nothing (MU.fromRows $ take i s1 ++ zipWith f (drop i s1) s2 ++ drop (n1 - i) s2)
    | otherwise = PWM Nothing (MU.fromRows $ take (-i) s2 ++ zipWith f (drop (-i) s2) s1 ++ drop (n2 + i) s1)
  where
    f = U.zipWith (\x y -> (x+y)/2)
    s1 = MU.toRows . _mat $ m1
    s2 = MU.toRows . _mat $ m2
    n1 = length s1
    n2 = length s2

mergePWMWeighted :: (PWM, [Int])  -- ^ pwm and weights at each position
                 -> (PWM, [Int])
                 -> Int                  -- ^ shift
                 -> (PWM, [Int])
mergePWMWeighted (pwm1, w1) (pwm2, w2) shift
    | shift >= 0 = mkPWM $ take shift rows1 ++
                   zipWith f (drop shift rows1) rows2 ++ drop (n1 - shift) rows2
    | otherwise = mkPWM $ take (-shift) rows2 ++
                  zipWith f (drop (-shift) rows2) rows1 ++ drop (n2 + shift) rows1
  where
    mkPWM xs = let (rs, ws) = unzip xs
               in (PWM Nothing $ MU.fromRows rs, ws)
    f (xs, a) (ys, b) = (U.zipWith (\x y -> (fromIntegral a * x + fromIntegral b * y) / fromIntegral (a + b)) xs ys, a + b)
    rows1 = zip (MU.toRows $ _mat pwm1) w1
    rows2 = zip (MU.toRows $ _mat pwm2) w2
    n1 = MU.rows $ _mat pwm1
    n2 = MU.rows $ _mat pwm2

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
            MS.SymMatrix _ vec <- MS.freeze mat
            let ((i,j), (d, (isSame, pos))) = V.minimumBy
                    (comparing (fst . snd)) $ V.map fromJust $ V.filter isJust vec
            if d < th
                then do
                    Just (nm1, pwm1, w1) <- VM.unsafeRead motifs' i
                    Just (nm2, pwm2, w2) <- VM.unsafeRead motifs' j
                    let merged = (nm1 ++ nm2, pwm', w')
                        (pwm',w') | isSame = mergePWMWeighted (pwm1, w1) (pwm2, w2) pos
                                  | otherwise = mergePWMWeighted (pwm1, w1)
                                              (rcPWM $ pwm2, reverse w2) pos

                    -- update
                    forM_ [0..n-1] $ \j' -> when (i /= j') $ do
                        x <- VM.unsafeRead motifs' j'
                        case x of
                            Just (_, pwm2',_) -> MSU.unsafeWrite mat (i,j') $ Just $
                                ((i,j'), align pwm' $ pwm2')
                            _ -> return ()
                    forM_ [0..n-1] $ \i' -> MSU.unsafeWrite mat (i',j) Nothing
                    VM.unsafeWrite motifs' i $ Just merged
                    VM.unsafeWrite motifs' j Nothing
                    iter mat
                else return ()

    -- initialization
    mat <- MSU.replicate (n,n) Nothing
    forM_ [0..n-1] $ \i -> forM_ [i+1 .. n-1] $ \j -> do
        Just (_, pwm1, _) <- VM.unsafeRead motifs' i
        Just (_, pwm2, _) <- VM.unsafeRead motifs' j
        MSU.unsafeWrite mat (i,j) $ Just $ ((i,j), align pwm1 pwm2)

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
