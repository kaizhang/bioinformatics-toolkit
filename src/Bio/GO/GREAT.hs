{-# LANGUAGE BangPatterns #-}
module Bio.GO.GREAT
    ( AssocRule(..)
    , getRegulatoryDomains
    , enrichedTerms
    ) where

import qualified Data.ByteString.Char8 as B
import Data.Conduit
import qualified Data.HashMap.Strict as M
import qualified Data.IntervalMap as IM
import Data.IntervalMap.Interval (lowerBound, upperBound)
import Data.Maybe (fromJust)
import Statistics.Distribution
import Statistics.Distribution.Hypergeometric (hypergeometric)

import Bio.Data.Bed
import Bio.GO

-- | how to associate genomic regions with genes
data AssocRule = BasalPlusExtension Int Int Int
               | TwoNearest
               | OneNearest

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
    go !n !m = do
        x <- await
        case x of
            Just bed ->do
                let chr = chrom bed
                    s = chromStart bed
                    e = chromEnd bed
                    terms = concatMap snd . IM.intersecting
                        (M.lookupDefault IM.empty chr tree) $ IM.IntervalCO s e
                go (n+1) $ foldl (\acc t -> M.insertWith (+) t 1 acc) m terms
            _ -> return (n, m)
{-# INLINE countTerms #-}

enrichedTerms :: (BEDLike b1, BEDLike b2, Monad m)
              => Source m b1
              -> Source m b2
              -> BEDTree [GOId]
              -> m (M.HashMap GOId Double)
enrichedTerms fg bg tree = do
    (n, c1) <- fg $$ countTerms tree
    (m, c2) <- bg $$ countTerms tree
    return $ flip M.mapWithKey c1 $ \t c ->
        let k = fromJust $ M.lookup t c2
        in complCumulative (hypergeometric k m n) $ fromIntegral c
