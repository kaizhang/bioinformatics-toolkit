{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE OverloadedStrings     #-}

module Bio.Data.Bed
    ( BEDLike(..)
    , BEDConvert(..)
    , BED
    , BED3
    , NarrowPeak
    , npSignal
    , npPvalue
    , npQvalue
    , npPeak
    , BEDExt(..)

    , BEDTree
    , bedToTree
    , sortedBedToTree
    , intersecting
    , isIntersected
    , sizeOverlapped
    , splitBed
    , splitBedBySize
    , splitBedBySizeLeft
    , splitBedBySizeOverlap
    , Sorted(..)
    , sortBed
    , intersectBed
    , intersectBedWith
    , intersectSortedBed
    , intersectSortedBedWith
    , isOverlapped
    , mergeBed
    , mergeBedWith
    , mergeSortedBed
    , mergeSortedBedWith
    , splitOverlapped
    , hReadBed
    , hReadBed'
    , readBed
    , readBed'
    , hWriteBed
    , hWriteBed'
    , writeBed
    , writeBed'

    -- * Utilities
    , fetchSeq
    , fetchSeq'
    , motifScan
    , getMotifScore
    , getMotifPValue
    , compareBed
    ) where

import           Conduit
import           Control.Arrow                ((***))
import           Control.Lens
import           Control.Monad.State.Strict
import qualified Data.ByteString.Char8        as B
import qualified Data.Foldable                as F
import           Data.Function                (on)
import qualified Data.HashMap.Strict          as M
import qualified Data.IntervalMap.Strict      as IM
import           Data.List                    (groupBy, sortBy)
import           Data.Maybe                   (fromJust)
import           Data.Ord                     (comparing)
import qualified Data.Vector                  as V
import qualified Data.Vector.Algorithms.Intro as I
import           System.IO

import           Bio.Data.Bed.Types
import Bio.Motif (Motif(..), Bkgd(..))
import qualified Bio.Motif                    as Motif
import qualified Bio.Motif.Search             as Motif
import           Bio.Seq
import           Bio.Seq.IO
import           Bio.Utils.Misc               (binBySize, binBySizeLeft,
                                               binBySizeOverlap, bins)

-- | Convert a set of sorted bed records to interval tree, with combining
-- function for equal keys.
sortedBedToTree :: (BEDLike b, F.Foldable f)
                => (a -> a -> a)
                -> Sorted (f (b, a))
                -> BEDTree a
sortedBedToTree f (Sorted xs) = M.fromList $
    map ((head *** IM.fromAscListWith f) . unzip) $ groupBy ((==) `on` fst) $
    map (\(b, x) -> (b^.chrom, (IM.IntervalCO (b^.chromStart) (b^.chromEnd), x))) $
    F.toList xs
{-# INLINE sortedBedToTree #-}

bedToTree :: BEDLike b
          => (a -> a -> a)
          -> [(b, a)]
          -> BEDTree a
bedToTree f xs = M.fromList $ map ((head *** IM.fromAscListWith f) . unzip) $
    groupBy ((==) `on` fst) $
    map (\(b, x) -> (b^.chrom, (IM.IntervalCO (b^.chromStart) (b^.chromEnd), x))) $
    V.toList $ V.create $ do
        v <- V.unsafeThaw $ V.fromList xs
        I.sortBy (compareBed `on` fst) v
        return v
{-# INLINE bedToTree #-}

intersecting :: BEDLike b => BEDTree a -> b -> IM.IntervalMap Int a
intersecting tree x = IM.intersecting (M.lookupDefault IM.empty (x^.chrom) tree) $
    IM.IntervalCO (x^.chromStart) $ x^.chromEnd
{-# INLINE intersecting #-}

isIntersected :: BEDLike b => BEDTree a -> b -> Bool
isIntersected tree = not . IM.null . intersecting tree
{-# INLINE isIntersected #-}

sizeOverlapped :: (BEDLike b1, BEDLike b2) => b1 -> b2 -> Int
sizeOverlapped b1 b2 | b1^.chrom /= b2^.chrom = 0
                     | overlap < 0 = 0
                     | otherwise = overlap
  where
    overlap = minimum [ b1^.chromEnd - b2^.chromStart
                      , b2^.chromEnd - b1^.chromStart
                      , b1^.chromEnd - b1^.chromStart
                      , b2^.chromEnd - b2^.chromStart ]

-- | split a bed region into k consecutive subregions, discarding leftovers
splitBed :: BEDConvert b => Int -> b -> [b]
splitBed k bed = map (uncurry (asBed (bed^.chrom))) $
    bins k (bed^.chromStart, bed^.chromEnd)
{-# INLINE splitBed #-}

-- | split a bed region into consecutive fixed size subregions, discarding leftovers
splitBedBySize :: BEDConvert b => Int -> b -> [b]
splitBedBySize k bed = map (uncurry (asBed (bed^.chrom))) $
    binBySize k (bed^.chromStart, bed^.chromEnd)
{-# INLINE splitBedBySize #-}

-- | split a bed region into consecutive fixed size subregions, including leftovers
splitBedBySizeLeft :: BEDConvert b => Int -> b -> [b]
splitBedBySizeLeft k bed = map (uncurry (asBed (bed^.chrom))) $
    binBySizeLeft k (bed^.chromStart, bed^.chromEnd)
{-# INLINE splitBedBySizeLeft #-}

splitBedBySizeOverlap :: BEDConvert b
                      => Int     -- ^ bin size
                      -> Int     -- ^ overlap size
                      -> b -> [b]
splitBedBySizeOverlap k o bed = map (uncurry (asBed (bed^.chrom))) $
    binBySizeOverlap k o (bed^.chromStart, bed^.chromEnd)
{-# INLINE splitBedBySizeOverlap #-}

-- | a type to imply that underlying data structure is sorted
newtype Sorted b = Sorted {fromSorted :: b} deriving (Show, Read, Eq)

-- | Compare bed records using only the chromosome, start and end positions.
-- Unlike the ``compare'' from the Ord type class, this function can compare
-- different types of BED data types.
compareBed :: (BEDLike b1, BEDLike b2) => b1 -> b2 -> Ordering
compareBed b1 b2 = compare (b1^.chrom, b1^.chromStart, b1^.chromEnd)
                           (b2^.chrom, b2^.chromStart, b2^.chromEnd)
{-# INLINE compareBed #-}

-- | sort BED, first by chromosome (alphabetical order), then by chromStart, last by chromEnd
sortBed :: BEDLike b => [b] -> Sorted (V.Vector b)
sortBed beds = Sorted $ V.create $ do
    v <- V.unsafeThaw . V.fromList $ beds
    I.sortBy compareBed v
    return v
{-# INLINE sortBed #-}

-- | return records in A that are overlapped with records in B
intersectBed :: (BEDLike b1, BEDLike b2, Monad m) => [b2] -> ConduitT b1 b1 m ()
intersectBed b = intersectSortedBed b'
  where
    b' = sortBed b
{-# INLINE intersectBed #-}

-- | return records in A that are overlapped with records in B
intersectSortedBed :: (BEDLike b1, BEDLike b2, Monad m)
                   => Sorted (V.Vector b2) -> ConduitT b1 b1 m ()
intersectSortedBed (Sorted b) = filterC (not . IM.null . intersecting tree)
  where
    tree = sortedBedToTree (\_ _ -> ()) . Sorted $ V.map (\x -> (x,())) b
{-# INLINE intersectSortedBed #-}

intersectBedWith :: (BEDLike b1, BEDLike b2, Monad m)
                 => ([b2] -> a)
                 -> [b2]
                 -> ConduitT b1 (b1, a) m ()
intersectBedWith fn = intersectSortedBedWith fn . sortBed
{-# INLINE intersectBedWith #-}

intersectSortedBedWith :: (BEDLike b1, BEDLike b2, Monad m)
                       => ([b2] -> a)
                       -> Sorted (V.Vector b2)
                       -> ConduitT b1 (b1, a) m ()
intersectSortedBedWith fn (Sorted b) = mapC f
  where
    f bed = (bed, fn $ concat $ IM.elems $ intersecting tree bed)
    tree = sortedBedToTree (++) $ Sorted $ V.map (\x -> (x, [x])) b
{-# INLINE intersectSortedBedWith #-}

isOverlapped :: (BEDLike b1, BEDLike b2) => b1 -> b2 -> Bool
isOverlapped b1 b2 = b1^.chrom == b2^.chrom &&
    not (b1^.chromEnd <= b2^.chromStart || b2^.chromEnd <= b1^.chromStart)

mergeBed :: (BEDConvert b, Monad m) => [b] -> ConduitT i b m ()
mergeBed = mergeSortedBed . sortBed
{-# INLINE mergeBed #-}

mergeBedWith :: (BEDLike b, Monad m)
             => ([b] -> a) -> [b] -> ConduitT i a m ()
mergeBedWith f = mergeSortedBedWith f . sortBed
{-# INLINE mergeBedWith #-}

mergeSortedBed :: (BEDConvert b, Monad m)
               => Sorted (V.Vector b)
               -> ConduitT i b m ()
mergeSortedBed = mergeSortedBedWith f
  where
    f xs = asBed (head xs ^. chrom) lo hi
      where
        lo = minimum $ map (^.chromStart) xs
        hi = maximum $ map (^.chromEnd) xs
{-# INLINE mergeSortedBed #-}

mergeSortedBedWith :: (BEDLike b, Monad m)
                   => ([b] -> a) -> Sorted (V.Vector b) -> ConduitT i a m ()
mergeSortedBedWith mergeFn (Sorted beds)
    | V.null beds = return ()
    | otherwise = do
        (_, r) <- V.foldM' f acc0 . V.tail $ beds
        yield $ mergeFn r
  where
    x0 = V.head beds
    acc0 = ((x0^.chrom, x0^.chromStart, x0^.chromEnd), [x0])
    f ((chr,lo,hi), acc) bed
        | chr /= chr' || s' > hi = yield (mergeFn acc) >>
                                   return ((chr',s',e'), [bed])
        | e' > hi = return ((chr',lo,e'), bed:acc)
        | otherwise = return ((chr,lo,hi), bed:acc)
      where
        chr' = bed^.chrom
        s' = bed^.chromStart
        e' = bed^.chromEnd
{-# INLINE mergeSortedBedWith #-}

-- | Split overlapped regions into non-overlapped regions. The input must be overlapped.
-- This function is usually used with `mergeBedWith`.
splitOverlapped :: BEDLike b => ([b] -> a) -> [b] -> [(BED3, a)]
splitOverlapped fun xs = filter ((>0) . size . fst) $
    evalState (F.foldrM f [] $ init xs') x0
  where
    x0 = (\(a,b) -> (fromEither a, M.singleton (b^.chromStart, b^.chromEnd) b)) $ last xs'
    xs' = sortBy (comparing (fromEither . fst)) $ concatMap
        ( \x -> [(Left $ x^.chromStart, x), (Right $ x^.chromEnd, x)] ) xs
    f (i, x) acc = do
        (j, set) <- get
        let bed = (asBed chr (fromEither i) j, fun $ M.elems set)
            set' = case i of
                Left _  -> M.delete (x^.chromStart, x^.chromEnd) set
                Right _ -> M.insert (x^.chromStart, x^.chromEnd) x set
        put (fromEither i, set')
        return (bed:acc)
    fromEither (Left x)  = x
    fromEither (Right x) = x
    chr = head xs ^. chrom
{-# INLINE splitOverlapped #-}

-- | Read records from a bed file handler in a streaming fashion.
hReadBed :: (BEDConvert b, MonadIO m) => Handle -> ConduitT i b m ()
hReadBed h = do eof <- liftIO $ hIsEOF h
                unless eof $ do
                    line <- liftIO $ B.hGetLine h
                    yield $ fromLine line
                    hReadBed h
{-# INLINE hReadBed #-}

-- | Non-streaming version.
hReadBed' :: (BEDConvert b, MonadIO m) => Handle -> m [b]
hReadBed' h = runConduit $ hReadBed h .| sinkList
{-# INLINE hReadBed' #-}

-- | Read records from a bed file in a streaming fashion.
readBed :: (BEDConvert b, MonadIO m) => FilePath -> ConduitT i b m ()
readBed fl = do handle <- liftIO $ openFile fl ReadMode
                hReadBed handle
                liftIO $ hClose handle
{-# INLINE readBed #-}

-- | Non-streaming version.
readBed' :: (BEDConvert b, MonadIO m) => FilePath -> m [b]
readBed' fl = runConduit $ readBed fl .| sinkList
{-# INLINE readBed' #-}

hWriteBed :: (BEDConvert b, MonadIO m) => Handle -> ConduitT b o m ()
hWriteBed handle = do
    x <- await
    case x of
        Nothing -> return ()
        Just bed -> (liftIO . B.hPutStrLn handle . toLine) bed >> hWriteBed handle
{-# INLINE hWriteBed #-}

hWriteBed' :: (BEDConvert b, MonadIO m) => Handle -> [b] -> m ()
hWriteBed' handle beds = runConduit $ yieldMany beds .| hWriteBed handle
{-# INLINE hWriteBed' #-}

writeBed :: (BEDConvert b, MonadIO m) => FilePath -> ConduitT b o m ()
writeBed fl = do handle <- liftIO $ openFile fl WriteMode
                 hWriteBed handle
                 liftIO $ hClose handle
{-# INLINE writeBed #-}

writeBed' :: (BEDConvert b, MonadIO m) => FilePath -> [b] -> m ()
writeBed' fl beds = runConduit $ yieldMany beds .| writeBed fl
{-# INLINE writeBed' #-}

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
