{-# LANGUAGE BangPatterns #-}
module Bio.GO.GREAT
    ( AssocRule(..)
    , getRegulatoryDomains
    , enrichedTerms
    ) where

import Control.Monad.Primitive
import qualified Data.ByteString.Char8 as B
import Conduit
import Data.Default.Class
import qualified Data.HashMap.Strict as M
import qualified Data.IntervalMap as IM
import qualified Data.Vector as V
import Data.Vector.Algorithms.Intro (sortBy)
import Data.List (sort, group)
import Data.Ord (comparing)
import Data.IntervalMap.Interval (lowerBound, upperBound)
import Data.Maybe

import Bio.Data.Bed
import Bio.GO
import Bio.Utils.Functions

-- | how to associate genomic regions with genes
data AssocRule = BasalPlusExtension Int Int Int
               | TwoNearest
               | OneNearest

instance Default AssocRule where
    def = BasalPlusExtension 5000 1000 50000

type Gene = ( B.ByteString  -- ^ chromosome
            , Int           -- ^ tss
            , Bool          -- ^ is forward stranded
            , [GOId]
            )

-- | given a gene list and the rule, compute the rulatory domain for each gene
getRegulatoryDomains :: AssocRule -> [Gene] -> BEDTree [GOId]
getRegulatoryDomains (BasalPlusExtension up dw ext) gs =
    bedToTree undefined $ flip map basal $ \(BED3 chr s e, go) ->
        let intervals = fromJust . M.lookup chr $ basalTree
            s' = let im = IM.intersecting intervals $ IM.OpenInterval (s-ext) s
                 in if IM.null im
                     then s - ext
                     else min s $ maximum $ map upperBound $ IM.keys im
            e' = let im = IM.intersecting intervals $ IM.OpenInterval e (e+ext)
                 in if IM.null im
                     then e + ext
                     else max e $ minimum $ map lowerBound $ IM.keys im
        in (BED3 chr s' e', go)
  where
    basalTree = bedToTree (error "encounter redundant regions") basal
    basal = flip map gs $ \(chr,tss,str,go) ->
        if str then (BED3 chr (tss - up) (tss + dw), go)
               else (BED3 chr (tss - dw) (tss + up), go)
getRegulatoryDomains _ _ = undefined
{-# INLINE getRegulatoryDomains #-}

-- | how many times a particular GO term is hit by given regions
countTerms :: (BEDLike b, Monad m)
           => BEDTree [GOId]
           -> Sink b m (Int, M.HashMap GOId Int)
countTerms tree = go 0 M.empty
  where
    go !n !termCount = do
        x <- await
        case x of
            Just bed ->do
                let chr = chrom bed
                    s = chromStart bed
                    e = chromEnd bed
                    terms = nub' . concat . IM.elems . IM.intersecting
                        (M.lookupDefault IM.empty chr tree) $ IM.IntervalCO s e
                go (n+1) $ foldl (\acc t -> M.insertWith (+) t 1 acc) termCount terms
            _ -> return (n, termCount)
{-
    expand = nub' . foldl f []
      where
        f acc i = traceBack i $ i : acc
        traceBack i acc = case getParentById i goMap of
            Just g -> traceBack (_oboId g) $ (_oboId g) : acc
            _ -> acc
-}
{-# INLINE countTerms #-}

nub' :: [B.ByteString] -> [B.ByteString]
nub' = map head . group . sort
{-# INLINE nub' #-}

{-
-- | since GO terms are organized as a tree, a node is considered as being hit
-- if any of its children is hit. During this step, we traverse GO tree to include
-- parents if the regions that hit parents are different from those hiting their
-- children
expandTermCounts :: M.HashMap GOId Int -> GOMap -> M.HashMap GOId Int
expandTermCounts counts goMap = M.fromList . prune . expand $ counts
  where
    prune m = loop goTree
      where
        loop (Node p children : xs) =
            let children' = filter (\(Node (_,c) _) -> c /= 0) children
                l = length children'
            in case () of
               _ | l == 0 -> p : loop xs
                 | snd p == (snd . rootLabel . head) children' -> loop $ children' ++ xs
                 | otherwise -> p : (loop $ children' ++ xs)
        loop _ = []

        goTree :: [Tree (GOId, Int)]
        goTree = map (fmap (\x -> (_oboId x, M.lookupDefault 0 (_oboId x) m))) . buildGOTree $ goMap
    expand m = M.foldlWithKey' f M.empty m
      where
        f acc i c = traceBack c i $ M.insertWith (+) i c acc
    traceBack c i m = case getParentById i goMap of
        Just g -> traceBack c (_oboId g) $ M.insertWith (+) (_oboId g) c m
        _ -> m
{-# INLINE expandTermCounts #-}
-}

enrichedTerms :: (BEDLike b1, BEDLike b2, PrimMonad m)
              => Source m b1
              -> Source m b2
              -> BEDTree [GOId]
              -> m (V.Vector (GOId, (Double, Double)))
enrichedTerms fg bg tree = do
    (n1, table1) <- fg $$ countTerms tree
    (n2, table2) <- bg $$ countTerms tree
    v <- V.unsafeThaw $ V.fromList $ M.toList $ flip M.mapWithKey table1 $ \t c ->
        let k = fromMaybe (error "x") $ M.lookup t table2
        in (1 - hyperquick c k n1 n2, fromIntegral (c*n2) / fromIntegral (n1*k))
    sortBy (comparing snd) v
    V.unsafeFreeze v
