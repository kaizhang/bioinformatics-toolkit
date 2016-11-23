{-# LANGUAGE BangPatterns #-}
module Bio.GO.GREAT
    ( AssocRule(..)
    , defaultAssocRule
    , getRegulatoryDomains
    , read3DContact
    , get3DRegulatoryDomains
    , enrichedTerms
    ) where

import           Conduit
import           Control.Monad                (forM_)
import           Control.Monad.Primitive
import qualified Data.ByteString.Char8        as B
import           Data.Function                (on)
import qualified Data.HashMap.Strict          as M
import qualified Data.IntervalMap             as IM
import           Data.List                    (foldl', group, sort, sortBy)
import           Data.Maybe
import           Data.Ord                     (comparing)
import qualified Data.Vector                  as V
import qualified Data.Vector.Algorithms.Intro as I

import           Bio.Data.Bed
import           Bio.GO
import           Bio.Utils.Functions
import           Bio.Utils.Misc               (readInt)

-- | How to associate genomic regions with genes
data AssocRule =
    BasalPlusExtension Int Int Int  -- ^ Upstream, downstream and extension
  | TwoNearest                      -- ^ Not implemented
  | OneNearest                      -- ^ Not implemented

-- | The default rule is defined as -5k ~ +1k around TSSs,
-- plus up to 1000k extension.
defaultAssocRule :: AssocRule
defaultAssocRule = BasalPlusExtension 5000 1000 1000000

-- | A Gene consists of the chromosome name, TSS, strandness and an associated value.
type Gene a = ((B.ByteString, Int, Bool), a)

-- | Given a gene list and the rule, compute the rulatory domain for each gene
getRegulatoryDomains :: AssocRule -> [Gene a] -> [(BED3, a)]
getRegulatoryDomains (BasalPlusExtension up dw ext) genes = (extendTail r ext, a) : rs
  where
    (rs, Just (r,a)) = foldl' f ([], Nothing) $ sortBy (compareBed `on` fst) basal
    f (acc, Nothing) (b,x) = (acc, Just (extendHead b ext, x))
    f (acc, Just (b', x')) (b,x)
        | chrom b' /= chrom b = ( (extendTail b' ext, x') : acc
                                , Just (extendHead b ext, x) )
        | chromEnd b' >= chromStart b = ((b',x') : acc, Just (b,x))
        | otherwise = let ext' = min ext $ (chromStart b - chromEnd b') `div` 2
                      in ((extendTail b' ext', x') : acc, Just (extendHead b ext', x))
    extendHead (BED3 chr s e) l | s - l >= 0 = BED3 chr (s-l) e
                                | otherwise = BED3 chr 0 e
    extendTail (BED3 chr s e) l = BED3 chr s (e+l)
    basal = flip map genes $ \((chr, tss, str), x) ->
        if str then (BED3 chr (tss - up) (tss + dw), x)
               else (BED3 chr (tss - dw) (tss + up), x)
getRegulatoryDomains _ _ = undefined
{-# INLINE getRegulatoryDomains #-}

-- | Read 3D contacts from a file, where each line contains 6 fields separated
-- by Tabs corresponding to 2 interacting loci. Example:
-- chr1 [TAB] 11 [TAB] 100 [TAB] chr2 [TAB] 23 [TAB] 200
read3DContact :: FilePath -> Source (ResourceT IO) (BED3, BED3)
read3DContact input = sourceFileBS input =$= linesUnboundedAsciiC =$= mapC f
  where
    f x = let (chr1:s1:e1:chr2:s2:e2:_) = B.split '\t' x
          in (BED3 chr1 (readInt s1) (readInt e1), BED3 chr2 (readInt s2) (readInt e2))

-- | 3D regulatory domains of a gene include the gene's promoter and regions
-- that contact with the gene's promoter.
get3DRegulatoryDomains :: Monad m
                       => [Gene a]   -- ^ Genes
                       -> Int        -- ^ Upstream
                       -> Int        -- ^ Downstream
                       -> Conduit (BED3, BED3) m (BED3, a)
get3DRegulatoryDomains genes up dw = concatRegions >> promoters
  where
    promoters = forM_ genes $ \((chr, tss, str), x) ->
        if str then yield (BED3 chr (gateZero $ tss - up) (tss + dw), x)
               else yield (BED3 chr (gateZero $ tss - dw) (tss + up), x)
    concatRegions = concatMapC $ \(locA, locB) ->
        map ((,) locB) (IM.elems $ intersecting basal locA) ++
        map ((,) locA) (IM.elems $ intersecting basal locB)
    basal = bedToTree undefined $ flip map genes $ \((chr, tss, str), x) ->
        if str then (BED3 chr (tss - up) (tss + dw), x)
               else (BED3 chr (tss - dw) (tss + up), x)
    gateZero x | x < 0 = 0
               | otherwise = x
{-# INLINE get3DRegulatoryDomains #-}

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
    I.sortBy (comparing snd) v
    V.unsafeFreeze v
