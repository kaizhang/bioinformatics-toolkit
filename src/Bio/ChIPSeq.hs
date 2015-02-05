{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
module Bio.ChIPSeq
    ( rpkmBed
    , rpkmSortedBed
    , profiling
    , rpkmBam
    , tagCountDistr
    ) where

import Bio.SamTools.Bam
import qualified Bio.SamTools.BamIndex as BI
import Control.Arrow ((***))
import Control.Monad (liftM, forM_)
import Control.Monad.Primitive (PrimMonad)
import Control.Monad.Trans.Class (lift)
import Data.Conduit
import qualified Data.Conduit.List as CL
import Data.Function (on)
import qualified Data.Foldable as F
import qualified Data.HashMap.Strict as M
import qualified Data.IntervalMap as IM
import Data.List (groupBy)
import Data.Maybe (fromJust, isJust)
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Algorithms.Intro as I
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.HashTable.IO as HT

import Bio.Data.Bam
import Bio.Data.Bed

-- | calculate RPKM on a set of unique regions. Regions (in bed format) would be kept in
-- memory but not tag file.
-- RPKM: Readcounts per kilobase per million reads. Only counts the starts of tags
rpkmBed :: (PrimMonad m, BEDLike b, G.Vector v Double)
     => [b] -> Sink BED m (v Double)
rpkmBed regions = do
    v <- lift $ do v' <- V.unsafeThaw . V.fromList . zip [0..] $ regions
                   I.sortBy (compareBed `on` snd) v'
                   V.unsafeFreeze v'
    let (idx, sortedRegions) = V.unzip v
        n = G.length idx
    rc <- rpkmSortedBed $ Sorted sortedRegions

    lift $ do
        result <- GM.new n
        G.sequence_ . G.imap (\x i -> GM.unsafeWrite result i (rc U.! x)) $ idx
        G.unsafeFreeze result
{-# INLINE rpkmBed #-}

-- | calculate RPKM on a set of regions. Regions must be sorted. The Sorted data
-- type is used to remind users to sort their data.
rpkmSortedBed :: (PrimMonad m, BEDLike b, G.Vector v Double)
              => Sorted (V.Vector b) -> Sink BED m (v Double)
rpkmSortedBed (Sorted regions) = do
    vec <- lift $ GM.replicate l 0
    n <- CL.foldM (f vec) (0 :: Int)
    let factor = fromIntegral n / 1e9
    lift $ liftM (G.imap (\i x -> x / factor / (fromIntegral . size) (regions V.! i)))
         $ G.unsafeFreeze vec
  where
    f v nTags tag = do
        let chr = chrom tag
            p | _strand tag == Just True = chromStart tag
              | _strand tag == Just False = chromEnd tag - 1
              | otherwise = error "Unkown strand"
            xs = snd . unzip $
                IM.containing (M.lookupDefault IM.empty chr intervalMap) p
        addOne v xs
        return $ succ nTags

    intervalMap = sortedBedToTree errMsg. Sorted . G.toList . G.zip regions . G.enumFromN 0 $ l
    addOne v' = mapM_ $ \x -> GM.unsafeRead v' x >>= GM.unsafeWrite v' x . (+1)
    l = G.length regions
    errMsg = error "rpkmSortedBed: redundant records"
{-# INLINE rpkmSortedBed #-}

-- | divide each region into consecutive bins, and count tags for each bin
profiling :: (PrimMonad m, G.Vector v Int, BEDLike b)
          => Int                   -- ^ bin size
          -> Sorted (V.Vector b)   -- ^ regions
          -> Sink BED m [v Int]
profiling k (Sorted beds) = do
    vectors <- lift $ G.forM beds $ \bed -> do
        let start = chromStart bed
            end = chromEnd bed
            num = (end - start) `div` k
            index i = (i - start) `div` k
        v <- GM.replicate num 0
        return (v, index)

    sink vectors
  where
    sink vs = do
        tag <- await
        case tag of
            Just (BED chr start end _ _ strand) -> do
                let p | strand == Just True = start
                      | strand == Just False = end - 1
                      | otherwise = error "unkown strand"
                    overlaps = snd . unzip $
                        IM.containing (M.lookupDefault IM.empty chr intervalMap) p
                lift $ forM_ overlaps $ \x -> do
                    let (v, f) = vs `G.unsafeIndex` x
                        i = f p
                    GM.unsafeRead v i >>= GM.unsafeWrite v i . (+1)
                sink vs

            _ -> lift $ mapM (G.unsafeFreeze . fst) $ G.toList vs
                                                            
    intervalMap = M.fromList
           . map ((head *** IM.fromAscListWith (error "profiling: non-unique regions")) . unzip)
           . groupBy ((==) `on` fst)
           . map (\(b, i) -> (chrom b, (IM.ClosedInterval (chromStart b) (chromEnd b), i)))
           . G.toList . G.zip beds
           $ G.enumFromN 0 n

    n = G.length beds
{-# INLINE profiling #-}

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
                           yield $ rc * 1e9 / n / fromIntegral (e-s)
                           conduit n h
    readCount l u = CL.foldM f 0.0
      where
        f acc bam = do let p1 = fromIntegral . fromJust . position $ bam
                           rl = fromIntegral . fromJust . queryLength $ bam
                           p2 = p1 + rl - 1
                       return $ if isReverse bam
                                   then if l <= p2 && p2 < u then acc + 1
                                                             else acc
                                   else if l <= p1 && p1 < u then acc + 1
                                                             else acc
{-# INLINE rpkmBam #-}

tagCountDistr :: G.Vector v Int => Sink BED IO (v Int)
tagCountDistr = loop M.empty
  where
    loop m = do
        x <- await
        case x of
            Just (BED chr s e _ _ (Just str)) -> do
                let p | str = s
                      | otherwise = 1 - e
                case M.lookup chr m of
                    Just table -> do
                        lift $ do c <- HT.lookup table p
                                  if isJust c 
                                     then HT.insert table p $ fromJust c + 1
                                     else HT.insert table p 1
                        loop m
                    _ -> do
                        t <- lift $ do t <- HT.new :: IO (HT.CuckooHashTable Int Int)
                                       HT.insert t p 1
                                       return t
                        loop $ M.insert chr t m
            _ -> lift $ do
                vec <- GM.replicate 100 0
                F.forM_ m $ \table ->
                    flip HT.mapM_ table $ \(_,v) -> do
                        let i = min 99 v
                        GM.unsafeRead vec i >>= GM.unsafeWrite vec i . (+1)
                G.unsafeFreeze vec
{-# INLINE tagCountDistr #-}

{-

tagDistrFromBamFile :: G.Vector v Int => FilePath -> IO (v Int)
tagDistrFromBamFile bamFl = do
    chrInfo <- readBam bamFl $= bamToBed $$ getRange
    readBam bamFl $= bamToBed $$ tagDistr chrInfo 

-- | the distribution of the number of tags per position
tagDistr :: G.Vector v Int => M.HashMap B.ByteString Int -> Sink BED IO (v Int)
tagDistr chrInfo = do
    [forwardCount,reverseCount] <- lift $ replicateM 2 $ forM chrInfo $ \v -> do
        ar <- newByteArray v
        fillByteArray ar 0 v 0
        return ar
    loop (forwardCount, reverseCount)
  where
    loop (f,r) = do
        x <- await
        case x of
            Just (BED chr s e _ _ (Just str)) ->
                if str
                   then do
                       let vec = M.lookupDefault (error "unknown chromosome") chr f
                       lift $ readByteArray vec s >>= writeByteArray vec s . (+(1::Word8))
                       loop (f,r)
                   else do
                       let vec = M.lookupDefault (error "unknown chromosome") chr r
                       lift $ readByteArray vec (e-1) >>= writeByteArray vec (e-1) . (+(1::Word8))
                       loop (f,r)
            _ -> lift $ do
                vec <- GM.replicate 100 0
                g f vec >> g r vec >> G.unsafeFreeze vec

    g xs v = F.forM_ xs $ \ar -> do
        let n = sizeofMutableByteArray ar
        forM_ [0..n-1] $ \i -> do
            c <- fmap (min 99 . fromIntegral) (readByteArray ar i :: IO Word8)
            GM.read v c >>= GM.write v c . (+1)
{-# INLINE tagDistr #-}

-- | loop over bed file and ranges for each chromosome
getRange :: Monad m => Sink BED m (M.HashMap B.ByteString Int)
getRange = go M.empty
  where
    go m = do
        x <- await
        case x of
            Just (BED chr _ e _ _ _) -> go $ M.insertWith max chr e m
            _ -> return m
{-# INLINE getRange #-}
-}
