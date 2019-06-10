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
    , _bed
    , _data

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

    -- * IO
    , streamBed
    , streamBedGzip
    , readBed
    , sinkFileBed
    , sinkFileBedGzip
    , sinkHandleBed
    , writeBed

    , compareBed
    ) where

import           Conduit
import           Control.Arrow                ((***))
import           Lens.Micro
import           Control.Monad.State.Strict
import qualified Data.ByteString.Char8        as B
import qualified Data.Foldable                as F
import           Data.Function                (on)
import qualified Data.HashMap.Strict          as M
import qualified Data.IntervalMap.Strict      as IM
import           Data.List                    (groupBy, sortBy)
import           Data.Ord                     (comparing)
import           Data.Conduit.Zlib           (gzip, ungzip, multiple)
import qualified Data.Vector                  as V
import qualified Data.Vector.Algorithms.Intro as I
import           System.IO

import           Bio.Data.Bed.Types
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
                 => (b1 -> [b2] -> a)
                 -> [b2]
                 -> ConduitT b1 a m ()
intersectBedWith fn = intersectSortedBedWith fn . sortBed
{-# INLINE intersectBedWith #-}

intersectSortedBedWith :: (BEDLike b1, BEDLike b2, Monad m)
                       => (b1 -> [b2] -> a)
                       -> Sorted (V.Vector b2)
                       -> ConduitT b1 a m ()
intersectSortedBedWith fn (Sorted b) = mapC $ \input -> fn input
    $ concat $ IM.elems $ intersecting tree input
  where
    tree = sortedBedToTree (++) $ Sorted $ V.map (\x -> (x, [x])) b
{-# INLINE intersectSortedBedWith #-}

isOverlapped :: (BEDLike b1, BEDLike b2) => b1 -> b2 -> Bool
isOverlapped b1 b2 = b1^.chrom == b2^.chrom &&
    not (b1^.chromEnd <= b2^.chromStart || b2^.chromEnd <= b1^.chromStart)

-- | Merge overlapping regions.
mergeBed :: (BEDConvert b, Monad m) => [b] -> ConduitT i b m ()
mergeBed xs = yieldMany xs' .| mergeSortedBed
  where
    Sorted xs' = sortBed xs
{-# INLINE mergeBed #-}

-- | Merge overlapping regions according to a merging function.
mergeBedWith :: (BEDLike b, Monad m)
             => ([b] -> a) -> [b] -> ConduitT i a m ()
mergeBedWith f xs = yieldMany xs' .| mergeSortedBedWith f
  where
    Sorted xs' = sortBed xs
{-# INLINE mergeBedWith #-}

-- | Merge overlapping regions. The input stream must be sorted first.
mergeSortedBed :: (BEDConvert b, Monad m) => ConduitT b b m ()
mergeSortedBed = mergeSortedBedWith f
  where
    f xs = asBed (head xs ^. chrom) lo hi
      where
        lo = minimum $ map (^.chromStart) xs
        hi = maximum $ map (^.chromEnd) xs
{-# INLINE mergeSortedBed #-}

-- | Merge overlapping regions according to a merging function. The input
-- stream must be sorted first.
mergeSortedBedWith :: (BEDLike b, Monad m)
                   => ([b] -> a) -> ConduitT b a m ()
mergeSortedBedWith mergeFn = headC >>= ( maybe mempty $ \b0 ->
    go ((b0^.chrom, b0^.chromStart, b0^.chromEnd), [b0]) )
  where
    go ((chr, s, e), acc) = headC >>= maybe (yield $ mergeFn acc) f
      where
        f bed | chr /= chr' || s' > e =
                    yield (mergeFn acc) >> go ((chr',s',e'), [bed])
              | s' < s = error "input stream is not sorted"
              | e' > e = go ((chr',s,e'), bed:acc)
              | otherwise = go ((chr,s,e), bed:acc)
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


streamBed :: (MonadResource m, BEDConvert b, MonadIO m)
          => FilePath -> ConduitT i b m () 
streamBed input = sourceFile input .| bsToBed
{-# INLINE streamBed #-}

streamBedGzip :: (BEDConvert b, MonadResource m, MonadThrow m, PrimMonad m)
              => FilePath -> ConduitT i b m () 
streamBedGzip input = sourceFile input .| multiple ungzip .| bsToBed
{-# INLINE streamBedGzip #-}

readBed :: BEDConvert b => FilePath -> IO [b]
readBed fl = runResourceT $ runConduit $ streamBed fl .| sinkList
{-# INLINE readBed #-}

sinkFileBed :: (BEDConvert b, MonadResource m) => FilePath -> ConduitT b o m ()
sinkFileBed output = bedToBS .| sinkFile output
{-# INLINE sinkFileBed #-}

sinkFileBedGzip :: (BEDConvert b, MonadResource m, MonadThrow m, PrimMonad m)
                => FilePath -> ConduitT b o m ()
sinkFileBedGzip output = bedToBS .| gzip .| sinkFile output
{-# INLINE sinkFileBedGzip #-}

sinkHandleBed :: (BEDConvert b, MonadIO m) => Handle -> ConduitT b o m ()
sinkHandleBed hdl = bedToBS .| sinkHandle hdl
{-# INLINE sinkHandleBed #-}

writeBed :: BEDConvert b => FilePath -> [b] -> IO ()
writeBed fl beds = runResourceT $ runConduit $ yieldMany beds .| sinkFileBed fl
{-# INLINE writeBed #-}

bedToBS :: (BEDConvert b, Monad m) => ConduitT b B.ByteString m ()
bedToBS = mapC toLine .| unlinesAsciiC
{-# INLINE bedToBS #-}

bsToBed :: (BEDConvert b, Monad m) => ConduitT B.ByteString b m ()
bsToBed = linesUnboundedAsciiC .| mapC fromLine
{-# INLINE bsToBed #-}
