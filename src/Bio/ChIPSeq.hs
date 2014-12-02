{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
module Bio.ChIPSeq
    ( rpkmBed
    , rpkmSortedBed
    , rpkmBam
    ) where

import Bio.Data.Bam
import Bio.Data.Bed
import Bio.SamTools.Bam
import qualified Bio.SamTools.BamIndex as BI
import Control.Monad.Primitive
import Control.Monad.Trans.Class (lift)
import Data.Conduit
import qualified Data.Conduit.List as CL
import Data.Function (on)
import qualified Data.HashMap.Strict as M
import qualified Data.IntervalMap as IM
import Data.Maybe
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Algorithms.Intro as I
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM

-- | calculate RPKM on a set of regions. Regions (in bed format) would be kept in
-- memory but not tag file.
-- RPKM: Readcounts per kilobase per million reads. Only counts the starts of tags
rpkmBed :: (PrimMonad m, BEDLike b, G.Vector v Double)
     => [b] -> Sink BED m (v Double)
rpkmBed regions = do
    v <- lift . V.unsafeThaw . V.fromList . zip [0..] $ regions
    lift $ I.sortBy (compareBed `on` snd) v
    v' <- lift $ V.unsafeFreeze v
    let (idx, sortedRegions) = V.unzip v'
        n = G.length idx
    rc <- rpkmSortedBed $ Sorted sortedRegions
    return $ G.create $ GM.new n >>= go n (rc::U.Vector Double) idx 0
  where
    go n' rc' idx' !i vec | i >= n' = return vec
                          | otherwise = do
                              let i' = idx' G.! i
                                  x = rc' G.! i
                              GM.unsafeWrite vec i' x
                              go n' rc' idx' (i+1) vec
{-# INLINE rpkmBed #-}

-- | calculate RPKM on a set of regions. Regions must be sorted. The Sorted data
-- type is used to remind users to sort their data.
rpkmSortedBed :: (PrimMonad m, BEDLike b, G.Vector v Double)
              => Sorted (V.Vector b) -> Sink BED m (v Double)
rpkmSortedBed (Sorted regions) = lift (GM.replicate n 0) >>= sink (0::Int)
  where
    sink !nTags v = do
        s <- await
        case s of
            Nothing -> lift $ do
                let f i bed = do
                        x <- GM.unsafeRead v i
                        GM.unsafeWrite v i $ x * 1e9
                                               / fromIntegral nTags
                                               / (fromIntegral . size) bed
                        return $ i + 1
                G.foldM'_ f 0 regions
                G.unsafeFreeze v
            Just tag -> do
                let chr = chrom tag
                    p | _strand tag == Just True = chromStart tag
                      | _strand tag == Just False = chromEnd tag
                      | otherwise = error "Unkown strand"
                    xs = snd . unzip $
                        IM.containing (M.lookupDefault IM.empty chr intervalMap) p
                lift $ addOne v xs
                sink (nTags+1) v
    intervalMap = M.fromList $ (c, IM.fromAscList x) : xs
      where
        (xs, (x, c, _)) = G.foldr f ([],([],"dummy",n-1)) regions
        f b (acc1, (acc2, !chr, !i))
            | chr == chr' = (acc1, (record : acc2, chr, i-1))
            | otherwise = ((chr, IM.fromAscList acc2) : acc1, ([record], chr', i-1))
          where
            record = (IM.ClosedInterval (chromStart b) (chromEnd b), i)
            chr' = chrom b
    addOne v' = mapM_ $ \x -> GM.unsafeRead v' x >>= GM.unsafeWrite v' x . (+1)
    n = G.length regions
{-# INLINE rpkmSortedBed #-}

-- | calculate RPKM using BAM file (*.bam) and its index file (*.bam.bai), using 
-- constant space
rpkmBam :: BEDLike b => FilePath -> Conduit b IO Double
rpkmBam fl = do
    nTags <- lift $ readBam fl $$ CL.foldM (\acc bam -> return $
                                  if isUnmap bam then acc else acc + 1) 0.0
    handle <- lift $ BI.open fl
    conduit nTags handle
  where
    conduit n h = do
        x <- await
        case x of
            Nothing -> lift $ BI.close h
            Just bed -> do let chr = chrom bed
                               s = chromStart bed
                               e = chromEnd bed
                           rc <- lift $ viewBam h (chr, s, e) $$ readCount s e
                           yield $ rc * 1e9 / n / fromIntegral (e-s+1)
                           conduit n h
    readCount l u = CL.foldM f 0.0
      where
        f acc bam = do let p1 = fromIntegral . fromJust . position $ bam
                           rl = fromIntegral . fromJust . queryLength $ bam
                           p2 = p1 + rl - 1
                       return $ if isReverse bam
                                   then if l <= p2 && p2 <= u then acc + 1
                                                              else acc
                                   else if l <= p1 && p1 <= u then acc + 1
                                                              else acc
{-# INLINE rpkmBam #-}
