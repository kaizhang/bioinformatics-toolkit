{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
module Bio.ChIPSeq
    ( rpkm
    , rpkm'
    ) where

import Bio.Data.Bed
import Control.Monad.Primitive
import Control.Monad.Trans.Class (lift)
import Data.Conduit
import qualified Data.HashMap.Strict as M
import qualified Data.IntervalMap as IM
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM

-- | calculate RPKM on a set of regions.
-- RPKM: Readcounts per kilobase per million reads. Only counts the starts of tags
rpkm :: (PrimMonad m, BEDFormat b, G.Vector v Double)
     => [b] -> Sink BED m (v Double)
rpkm regions = rpkm' sortedRegions
  where sortedRegions = sortBed . V.fromList $ regions
{-# INLINE rpkm #-}

-- | calculate RPKM on a set of regions. Regions must be sorted. The Sorted data
-- type is used to remind users to sort their data.
rpkm' :: (PrimMonad m, BEDFormat b, G.Vector v Double)
      => Sorted (V.Vector b) -> Sink BED m (v Double)
rpkm' (Sorted regions) = do
    vec <- lift $ GM.replicate n 0
    sink vec (0::Int)
  where
    sink v !nTags = do
        s <- await
        case s of
            Nothing -> lift $ do
                let f i bed = do
                        x <- GM.unsafeRead v i
                        GM.unsafeWrite v i $ x / ((fromIntegral . size) bed / 1000)
                                               / r
                        return $ i + 1
                    r = fromIntegral nTags / 1000000
                G.foldM'_ f 0 regions
                G.unsafeFreeze v
            Just tag -> do
                let chr = chrom tag
                    p | _strand tag == Just True = chromStart tag
                      | _strand tag == Just False = chromEnd tag
                      | otherwise = error "Unkown strand"
                    xs = snd . unzip $
                        IM.containing (M.lookupDefault errMsg chr intervalMap) p
                lift $ addOne v xs
                sink v (nTags+1)
    intervalMap = M.fromList . fst $ G.foldr f ([],([],"dummy",0)) regions
      where
        f b (acc1, (acc2, !chr, !i))
            | chr == chr' = (acc1, (record : acc2, chr, i+1))
            | otherwise = ((chr, IM.fromAscList acc2) : acc1, ([record], chr', i+1))
          where
            record = (IM.ClosedInterval (chromStart b) (chromEnd b), i)
            chr' = chrom b
    addOne v' = mapM_ $ \x -> GM.unsafeRead v' x >>= GM.unsafeWrite v' x . (+1)
    n = G.length regions
    errMsg = undefined
{-# INLINE rpkm' #-}
