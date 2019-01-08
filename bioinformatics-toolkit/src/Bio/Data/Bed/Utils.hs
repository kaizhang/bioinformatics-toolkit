{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE BangPatterns #-}

module Bio.Data.Bed.Utils
    ( fetchSeq
    , fetchSeq'
    , motifScan
    , getMotifScore
    , getMotifPValue
    , monoColonalize
    , BaseMap
    , baseMap
    , queryBaseMap
    , rpkmBed
    , rpkmSortedBed
    , countTagsBinBed
    , countTagsBinBed'
    , tagCountDistr
    , peakCluster
    ) where

import           Conduit
import           Control.Lens
import           Control.Monad.State.Strict
import qualified Data.ByteString.Char8 as B
import qualified Data.Foldable                as F
import           Data.Function                (on)
import qualified Data.HashMap.Strict          as M
import qualified Data.IntervalMap.Strict      as IM
import           Data.Maybe                   (fromJust, fromMaybe)
import qualified Data.Vector                  as V
import qualified Data.Vector.Algorithms.Intro as I
import qualified Data.Vector.Generic          as G
import qualified Data.Vector.Generic.Mutable  as GM
import qualified Data.Vector.Unboxed          as U

import           Bio.Data.Bed
import           Bio.Data.Bed.Types
import           Bio.Motif                    (Bkgd (..), Motif (..))
import qualified Bio.Motif                    as Motif
import qualified Bio.Motif.Search             as Motif
import           Bio.Seq hiding (length)
import           Bio.Seq.IO
import qualified Bio.Utils.BitVector as BV


-- | retreive sequences
fetchSeq :: (BioSeq DNA a, MonadIO m)
         => Genome
         -> ConduitT BED (Either String (DNA a)) m ()
fetchSeq g = mapMC f
  where
    f bed = do
        dna <- liftIO $ getSeq g (bed^.chrom, bed^.chromStart, bed^.chromEnd)
        return $ case bed^.strand of
            Just False -> rc <$> dna
            _          -> dna
{-# INLINE fetchSeq #-}

fetchSeq' :: (BioSeq DNA a, MonadIO m) => Genome -> [BED] -> m [Either String (DNA a)]
fetchSeq' g beds = runConduit $ yieldMany beds .| fetchSeq g .| sinkList
{-# INLINE fetchSeq' #-}

-- | Identify motif binding sites
motifScan :: (BEDLike b, MonadIO m)
          => Genome -> [Motif] -> Bkgd -> Double -> ConduitT b BED m ()
motifScan g motifs bg p = awaitForever $ \bed -> do
    r <- liftIO $ getSeq g (bed^.chrom, bed^.chromStart, bed^.chromEnd)
    case r of
        Left _    -> return ()
        Right dna -> mapM_ (getTFBS dna (bed^.chrom, bed^.chromStart)) motifs'
  where
    getTFBS dna (chr, s) (nm, (pwm, cutoff), (pwm', cutoff')) = toProducer
        ( (Motif.findTFBS bg pwm (dna :: DNA IUPAC) cutoff True .|
            mapC (\i -> bed & chromStart +~ i & chromEnd +~ i & strand .~ Just True)) >>
          (Motif.findTFBS bg pwm' dna cutoff' True .|
            mapC (\i -> bed & chromStart +~ i & chromEnd +~ i & strand .~ Just False)) )
      where
        n = Motif.size pwm
        bed = asBed chr s (s+n) & name .~ Just nm
    motifs' = flip map motifs $ \(Motif nm pwm) ->
        let cutoff = Motif.pValueToScore p bg pwm
            cutoff' = Motif.pValueToScore p bg pwm'
            pwm' = Motif.rcPWM pwm
        in (nm, (pwm, cutoff), (pwm', cutoff'))
{-# INLINE motifScan #-}

-- | Retrieve motif matching scores
getMotifScore :: MonadIO m
              => Genome -> [Motif] -> Bkgd -> ConduitT BED BED m ()
getMotifScore g motifs bg = awaitForever $ \bed -> do
    r <- liftIO $ getSeq g (bed^.chrom, bed^.chromStart, bed^.chromEnd)
    let r' = case bed^.strand of
            Just False -> rc <$> r
            _          -> r
    case r' of
        Left _ -> return ()
        Right dna -> do
            let pwm = M.lookupDefault (error "can't find motif with given name")
                        (fromJust $ bed^.name) motifMap
                sc = Motif.score bg pwm (dna :: DNA IUPAC)
            yield $ score .~ Just sc $ bed
  where
    motifMap = M.fromListWith (error "found motif with same name") $
        map (\(Motif nm pwm) -> (nm, pwm)) motifs
{-# INLINE getMotifScore #-}

getMotifPValue :: Monad m
               => Maybe Double   -- ^ whether to truncate the motif score CDF.
                                 -- Doing this will significantly reduce memory
                                 -- usage without sacrifice accuracy.
               -> [Motif] -> Bkgd -> ConduitT BED BED m ()
getMotifPValue truncation motifs bg = mapC $ \bed ->
    let nm = fromJust $ bed^.name
        sc = fromJust $ bed^.score
        d = M.lookupDefault (error "can't find motif with given name")
                nm motifMap
        p = 1 - Motif.cdf d sc
     in score .~ Just p $ bed
  where
    motifMap = M.fromListWith (error "getMotifPValue: found motif with same name") $
        map (\(Motif nm pwm) -> (nm, compressCDF $ Motif.scoreCDF bg pwm)) motifs
    compressCDF = case truncation of
        Nothing -> id
        Just x  -> Motif.truncateCDF x
{-# INLINE getMotifPValue #-}

-- | process a sorted BED stream, keep only mono-colonal tags
monoColonalize :: Monad m => ConduitT BED BED m ()
monoColonalize = do
    x <- headC
    case x of
        Just b -> yield b >> concatMapAccumC f b
        Nothing -> return ()
  where
    f cur prev = case compareBed prev cur of
        GT -> error $
            "Input is not sorted: " ++ show prev ++ " > " ++ show cur
        LT -> (cur, [cur])
        _ -> if prev^.strand == cur^.strand then (cur, []) else (cur, [cur])
{-# INLINE monoColonalize #-}

newtype BaseMap = BaseMap (M.HashMap B.ByteString BV.BitVector)

-- | Count the tags (starting positions) at each position in the genome.
baseMap :: PrimMonad m
        => [(B.ByteString, Int)]   -- ^ chromosomes and their sizes
        -> ConduitT BED o m BaseMap
baseMap chrs = do
    bvs <- lift $ fmap M.fromList $ forM chrs $ \(chr, n) -> do
        bv <- BV.zeros n
        return (chr, bv)

    mapM_C $ \bed -> case M.lookup (bed^.chrom) bvs of
        Nothing -> return ()
        Just bv -> if fromMaybe True $ bed^.strand
            then BV.set bv $ bed^.chromStart
            else BV.set bv $ bed^.chromEnd

    lift $ fmap BaseMap $ sequence $ fmap BV.unsafeFreeze bvs 

queryBaseMap :: BEDLike b => b -> BaseMap -> Maybe [Bool]
queryBaseMap bed (BaseMap bm) = case M.lookup (bed^.chrom) bm of
    Nothing -> Nothing
    Just bv ->
        let res = map (bv BV.!) [bed^.chromStart .. bed^.chromEnd - 1]
        in case bed^.strand of
            Just False -> Just $ reverse res
            _ -> Just res

-- | calculate RPKM on a set of unique regions. Regions (in bed format) would be kept in
-- memory but not tag file.
-- RPKM: Readcounts per kilobase per million reads. Only counts the starts of tags
rpkmBed :: (PrimMonad m, BEDLike b, G.Vector v Double)
     => [b] -> ConduitT BED o m (v Double)
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
              => Sorted (V.Vector b) -> ConduitT BED o m (v Double)
rpkmSortedBed (Sorted regions) = do
    vec <- lift $ GM.replicate l 0
    n <- foldMC (count vec) (0 :: Int)
    let factor = fromIntegral n / 1e9
    lift $ liftM (G.imap (\i x -> x / factor / (fromIntegral . size) (regions V.! i)))
         $ G.unsafeFreeze vec
  where
    count v nTags tag = do
        let p | tag^.strand == Just True = tag^.chromStart
              | tag^.strand == Just False = tag^.chromEnd - 1
              | otherwise = error "Unkown strand"
            xs = concat $ IM.elems $
                IM.containing (M.lookupDefault IM.empty (tag^.chrom) intervalMap) p
        addOne v xs
        return $ succ nTags

    intervalMap = sortedBedToTree (++) . Sorted . G.toList . G.zip regions .
                  G.map return . G.enumFromN 0 $ l
    addOne v' = mapM_ $ \x -> GM.unsafeRead v' x >>= GM.unsafeWrite v' x . (+1)
    l = G.length regions
{-# INLINE rpkmSortedBed #-}

-- | divide each region into consecutive bins, and count tags for each bin and
-- return the number of all tags. Note: a tag is considered to be overlapped
-- with a region only if the starting position of the tag is in the region. For
-- the common sense overlapping, use countTagsBinBed'.
countTagsBinBed :: (Integral a, PrimMonad m, G.Vector v a, BEDLike b)
           => Int   -- ^ bin size
           -> [b]   -- ^ regions
           -> ConduitT BED o m ([v a], Int)
countTagsBinBed k beds = do
    initRC <- lift $ forM beds $ \bed -> do
        let start = bed^.chromStart
            end = bed^.chromEnd
            num = (end - start) `div` k
            index i = (i - start) `div` k
        v <- GM.replicate num 0
        return (v, index)

    sink 0 $ V.fromList initRC
  where
    sink !nTags vs = do
        tag <- await
        case tag of
            Just bed -> do
                let p | bed^.strand == Just True = bed^.chromStart
                      | bed^.strand == Just False = bed^.chromEnd - 1
                      | otherwise = error "profiling: unkown strand"
                    overlaps = concat $ IM.elems $
                        IM.containing (M.lookupDefault IM.empty (bed^.chrom) intervalMap) p
                lift $ forM_ overlaps $ \x -> do
                    let (v, idxFn) = vs `G.unsafeIndex` x
                        i = let i' = idxFn p
                                l = GM.length v
                            in if i' >= l then l - 1 else i'
                    GM.unsafeRead v i >>= GM.unsafeWrite v i . (+1)
                sink (nTags+1) vs

            _ -> do rc <- lift $ mapM (G.unsafeFreeze . fst) $ G.toList vs
                    return (rc, nTags)

    intervalMap = bedToTree (++) $ zip beds $ map return [0..]
{-# INLINE countTagsBinBed #-}

-- | Same as countTagsBinBed, except that tags are treated as complete intervals
-- instead of single points.
countTagsBinBed' :: (Integral a, PrimMonad m, G.Vector v a, BEDLike b1, BEDLike b2)
                 => Int   -- ^ bin size
                 -> [b1]   -- ^ regions
                 -> ConduitT b2 o m ([v a], Int)
countTagsBinBed' k beds = do
    initRC <- lift $ forM beds $ \bed -> do
        let start = bed^.chromStart
            end = bed^.chromEnd
            num = (end - start) `div` k
            index i = (i - start) `div` k
        v <- GM.replicate num 0
        return (v, index)

    sink 0 $ V.fromList initRC
  where
    sink !nTags vs = do
        tag <- await
        case tag of
            Just bed -> do
                let chr = bed^.chrom
                    start = bed^.chromStart
                    end = bed^.chromEnd
                    overlaps = concat $ IM.elems $ IM.intersecting
                        (M.lookupDefault IM.empty chr intervalMap) $ IM.IntervalCO start end
                lift $ forM_ overlaps $ \x -> do
                    let (v, idxFn) = vs `G.unsafeIndex` x
                        lo = let i = idxFn start
                             in if i < 0 then 0 else i
                        hi = let i = idxFn end
                                 l = GM.length v
                             in if i >= l then l - 1 else i
                    forM_ [lo..hi] $ \i ->
                        GM.unsafeRead v i >>= GM.unsafeWrite v i . (+1)
                sink (nTags+1) vs

            _ -> do rc <- lift $ mapM (G.unsafeFreeze . fst) $ G.toList vs
                    return (rc, nTags)

    intervalMap = bedToTree (++) $ zip beds $ map return [0..]
{-# INLINE countTagsBinBed' #-}


{-
-- | calculate RPKM using BAM file (*.bam) and its index file (*.bam.bai), using
-- constant space
rpkmBam :: BEDLike b => FilePath -> Conduit b IO Double
rpkmBam fl = do
    nTags <- lift $ readBam fl $$ foldMC (\acc bam -> return $
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
    readCount l u = foldMC f 0.0
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
-}

tagCountDistr :: PrimMonad m => G.Vector v Int => ConduitT BED o m (v Int)
tagCountDistr = loop M.empty
  where
    loop m = do
        x <- await
        case x of
            Just bed -> do
                let p | fromMaybe True (bed^.strand) = bed^.chromStart
                      | otherwise = 1 - bed^.chromEnd
                case M.lookup (bed^.chrom) m of
                    Just table -> loop $ M.insert (bed^.chrom) (M.insertWith (+) p 1 table) m
                    _ -> loop $ M.insert (bed^.chrom) (M.fromList [(p,1)]) m
            _ -> lift $ do
                vec <- GM.replicate 100 0
                F.forM_ m $ \table ->
                    F.forM_ table $ \v -> do
                        let i = min 99 v
                        GM.unsafeRead vec i >>= GM.unsafeWrite vec i . (+1)
                G.unsafeFreeze vec
{-# INLINE tagCountDistr #-}

-- | cluster peaks
peakCluster :: (BEDLike b, Monad m)
            => [b]   -- ^ peaks
            -> Int   -- ^ radius
            -> Int   -- ^ cutoff
            -> ConduitT o BED m ()
peakCluster peaks r th = mergeBedWith mergeFn peaks' .| filterC g
  where
    peaks' = map f peaks
    f b = let c = (b^.chromStart + b^.chromEnd) `div` 2
          in asBed (b^.chrom) (c-r) (c+r) :: BED3
    mergeFn xs = asBed (head xs ^. chrom) lo hi & score .~ Just (fromIntegral $ length xs)
      where
        lo = minimum $ map (^.chromStart) xs
        hi = maximum $ map (^.chromEnd) xs
    g b = fromJust (b^.score) >= fromIntegral th
{-# INLINE peakCluster #-}
