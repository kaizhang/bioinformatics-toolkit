{-# LANGUAGE OverloadedStrings, UnicodeSyntax, BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}

module Bio.Utils.Transform (
      ihs
    , ihs'
    , copula
) where

import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import Statistics.Distribution
import Statistics.Distribution.Normal
import Statistics.Sample
import Data.List
import Control.Monad.ST
import Data.Function
import Control.Arrow

(!) ∷ G.Vector v a ⇒ v a → Int → a
(!) = (G.!)

-- | Inverse Hyperbolic sine transformation
ihs ∷ Double → Double → Double
ihs !θ !x | θ == 0 = x
          | otherwise = log (θ * x + sqrt (θ^(2∷Int) * x^(2∷Int) + 1)) / θ

ihs' ∷ Double → Double
ihs' = ihs 1

empiricalCDF ∷ G.Vector v Double ⇒ v Double → v Double
{-# INLINE empiricalCDF #-}
empiricalCDF xs = runST $ do
    let n = G.length xs
        indices = groupBy ( (==) `on` ((xs!).snd) ) $ zip [1.0..] $ sortBy (compare `on` (xs!)) [0..n-1]
        updates mv (v,ys) = mapM_ (flip (GM.unsafeWrite mv) v.snd) ys
    xs' ← G.thaw xs
    mapM_ (updates xs'. ((flip (/) (fromIntegral n).fst.last) &&& id)) indices
    G.unsafeFreeze xs'

cdf ∷ G.Vector v Double ⇒ v Double → v Double
{-# INLINE cdf #-}
cdf xs = let f = empiricalCDF xs
             n = fromIntegral $ G.length xs
             δ = 1 / (4 * (n**0.25) * sqrt (pi * log n))
         in G.map (\ x → case () of
             _ | x < δ → δ
               | x > 1 - δ → 1 - δ
               | otherwise → x) f

h ∷ G.Vector v Double ⇒ v Double → v Double
{-# INLINE h #-}
h = G.map (quantile standard) . cdf

copula ∷ G.Vector v Double ⇒ v Double → v Double
copula xs = let (μ, v) = meanVarianceUnb xs
                σ = sqrt v
            in G.map (\ x → μ + σ*x) $ h xs
