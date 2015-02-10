{-# LANGUAGE BangPatterns #-}
module Bio.GO.GREAT
    ( AssocRule(..)
    , getRegulatoryDomains
    , enrichedTerms
    ) where

import Control.Monad.Primitive
import qualified Data.ByteString.Char8 as B
import Data.Conduit
import Data.Default.Class
import qualified Data.HashMap.Strict as M
import qualified Data.IntervalMap as IM
import qualified Data.Vector as V
import Data.Vector.Algorithms.Intro (sortBy)
import Data.List (sort, group)
import Data.Ord (comparing)
import Data.IntervalMap.Interval (lowerBound, upperBound)
import Data.Maybe (fromJust)

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
            s' = case IM.intersecting intervals (IM.OpenInterval (s - ext)  s) of
                [] -> s - ext
                x -> min s (maximum $ map (upperBound . fst) x)
            e' = case IM.intersecting intervals (IM.OpenInterval e (e + ext)) of
                [] -> e + ext
                x -> max e (minimum $ map (lowerBound . fst) x)
        in (BED3 chr s' e', go)
  where
    basalTree = bedToTree (error "encounter redundant regions") basal
    basal = flip map gs $ \(chr,tss,str,go) ->
        if str then (BED3 chr (tss - up) (tss + dw), go)
               else (BED3 chr (tss - dw) (tss + up), go)
getRegulatoryDomains _ _ = undefined
{-# INLINE getRegulatoryDomains #-}

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
                    terms = nub' . concatMap snd . IM.intersecting
                        (M.lookupDefault IM.empty chr tree) $ IM.IntervalCO s e
                go (n+1) $ foldl (\acc t -> M.insertWith (+) t 1 acc) termCount terms
            _ -> return (n, termCount)
{-# INLINE countTerms #-}

nub' :: [B.ByteString] -> [B.ByteString]
nub' = map head . group . sort
{-# INLINE nub' #-}

enrichedTerms :: (BEDLike b1, BEDLike b2, PrimMonad m)
              => Source m b1
              -> Source m b2
              -> BEDTree [GOId]
              -> m (V.Vector (GOId, Double))
enrichedTerms fg bg tree = do
    (n1, table1) <- fg $$ countTerms tree
    (n2, table2) <- bg $$ countTerms tree
    v <- V.unsafeThaw $ V.fromList $ M.toList $ flip M.mapWithKey table1 $ \t c ->
        let k = fromJust $ M.lookup t table2
        in 1 - hyperquick c k n1 n2
    sortBy (comparing snd) v
    V.unsafeFreeze v
