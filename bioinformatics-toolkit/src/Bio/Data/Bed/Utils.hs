{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE BangPatterns #-}

module Bio.Data.Bed.Utils
    ( fetchSeq
    , clipBed
    , CutoffMotif(..)
    , mkCutoffMotif
    , scanMotif 
    , monoColonalize
    , BaseMap(..)
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
import           Lens.Micro
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
import System.IO

import           Bio.Data.Bed
import           Bio.Data.Bed.Types
import           Bio.Motif                    (Bkgd (..), Motif (..), CDF, PWM)
import qualified Bio.Motif                    as Motif
import qualified Bio.Motif.Search             as Motif
import           Bio.Seq hiding (length)
import           Bio.Seq.IO
import qualified Bio.Utils.BitVector as BV

clipBed :: (BEDLike b, Monad m)
        => [(B.ByteString, Int)]  -- ^ Chromosome sizes
        -> ConduitT b b m ()
clipBed chrsize = mapC f
  where
    f x = case M.lookup (x^.chrom) chrsize' of
        Nothing -> x
        Just n -> chromStart %~ max 0 $ chromEnd %~ min n $ x
    chrsize' = M.fromListWith (error "redundant chromosomes") chrsize
{-# INLINE clipBed #-}

-- | retreive sequences
fetchSeq :: BioSeq DNA a
         => Genome
         -> BED
         -> IO (Either String (DNA a))
fetchSeq g bed = do
    dna <- getSeq g (bed^.chrom, bed^.chromStart, bed^.chromEnd)
    return $ case bed^.strand of
        Just False -> rc <$> dna
        _          -> dna
{-# INLINE fetchSeq #-}

-- | Motif with predefined cutoff score. All necessary intermediate data
-- structure for motif scanning are stored.
data CutoffMotif = CutoffMotif
    { _motif_name :: B.ByteString
    , _motif_pwm :: PWM
    , _motif_sigma :: U.Vector Double
    , _motif_pwm_rc :: PWM
    , _motif_sigma_rc :: U.Vector Double
    , _background :: Bkgd
    , _cutoff :: Double
    , _cdf :: CDF }

mkCutoffMotif :: Bkgd
              -> Double  -- ^ p-value
              -> Motif -> CutoffMotif
mkCutoffMotif bg p motif = CutoffMotif (_name motif) (_pwm motif) sigma pwm'
    sigma' bg sc $ Motif.truncateCDF (1 - p * 10) cdf
  where
    cdf = Motif.scoreCDF bg $ _pwm motif
    sc = Motif.cdf' cdf $ 1 - p
    sigma = Motif.optimalScoresSuffix bg $ _pwm motif
    pwm' = Motif.rcPWM $ _pwm motif
    sigma' = Motif.optimalScoresSuffix bg pwm'

-- | Motif score is in [0, 1000]: ( 1 / (1 + exp (-(-logP - 5))) ) * 1000.
scanMotif :: (BEDLike b, MonadIO m)
          => Genome -> [CutoffMotif] -> ConduitT b BED m ()
scanMotif g motifs = awaitForever $ \bed -> do
    let (chr, start, end) = (bed^.chrom, bed^.chromStart, bed^.chromEnd)
    liftIO (getSeq g (chr, start, end)) >>= \case
        Left _    -> liftIO $ hPutStrLn stderr $
            "Warning: no sequence for region: " ++ show (chr, start, end)
        Right dna -> forM_ motifs $ \CutoffMotif{..} -> do
            let mkBed str (i, sc) = BED chr (start + i) (start + i + n)
                    (Just $ _motif_name) (Just $ toAffinity $ 1 - Motif.cdf _cdf sc)
                    (Just str)
                n = Motif.size _motif_pwm
            -- Scan forward strand
            Motif.findTFBSWith _motif_sigma _background _motif_pwm
                (dna :: DNA IUPAC) _cutoff True .| mapC (mkBed True)
            -- Scan reverse strand
            Motif.findTFBSWith _motif_sigma_rc _background _motif_pwm_rc
                dna _cutoff True .| mapC (mkBed False)
  where
    toAffinity x' = round $ sc * 1000
      where
        sc = 1 / (1 + exp (-(x - 5)))
        x = negate $ logBase 10 $ max 1e-20 x'
{-# INLINE scanMotif #-}

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
            else BV.set bv $ bed^.chromEnd - 1

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
    vs <- lift $ fmap V.fromList $ forM beds $ \bed -> do
        let start = bed^.chromStart
            num = ((bed^.chromEnd) - start) `div` k
            index i = (i - start) `div` k
        v <- GM.replicate num 0
        return (v, index)
    nTags <- foldMC (f vs) 0
    rc <- lift $ mapM (G.unsafeFreeze . fst) $ G.toList vs
    return (rc, nTags)
  where
    f vs n bed = do
        let pos | bed^.strand == Just True = bed^.chromStart
                | bed^.strand == Just False = bed^.chromEnd - 1
                | otherwise = error "unkown strand."
            overlaps = concat $ IM.elems $ IM.containing
                (M.lookupDefault IM.empty (bed^.chrom) intervalMap) pos
        forM_ overlaps $ \x -> do
            let (v, idxFn) = vs `G.unsafeIndex` x
                i = let i' = idxFn pos
                        l = GM.length v
                    in if i' >= l then l - 1 else i'
            GM.unsafeModify v (+1) i
        return $ n + 1
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
