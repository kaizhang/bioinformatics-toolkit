{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
--------------------------------------------------------------------------------
-- |
-- Module      :  $Header$
-- Copyright   :  (c) 2014 Kai Zhang
-- License     :  MIT

-- Maintainer  :  kai@kzhang.org
-- Stability   :  experimental
-- Portability :  portable

-- functions for processing BED files
--------------------------------------------------------------------------------

module Bio.Data.Bed
    ( BEDLike(..)
    , BEDTree
    , bedToTree
    , sortedBedToTree
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
    -- * BED6 format
    , BED(..)
    -- * BED3 format
    , BED3(..)
    , fetchSeq
    , fetchSeq'
    , compareBed
    , convert
    ) where

import Control.Arrow ((***))
import Control.Monad.State.Strict
import qualified Data.ByteString.Char8 as B
import Conduit
import Data.Default.Class (Default(..))
import Data.Function (on)
import qualified Data.Foldable as F
import qualified Data.HashMap.Strict as M
import qualified Data.IntervalMap.Strict as IM
import Data.List (groupBy)
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import qualified Data.Vector.Algorithms.Intro as I
import System.IO

import Bio.Seq
import Bio.Seq.IO
import Bio.Utils.Misc ( binBySizeLeft, binBySize, binBySizeOverlap, bins
                      , readInt, readDouble )

-- | a class representing BED-like data, e.g., BED3, BED6 and BED12. BED format
-- uses 0-based index (see documentation).
class BEDLike b where
    -- | construct bed record from chromsomoe, start location and end location
    asBed :: B.ByteString -> Int -> Int -> b

    -- | convert bytestring to bed format
    fromLine :: B.ByteString -> b

    -- | convert bed to bytestring
    toLine :: b -> B.ByteString

    -- | field accessor
    chrom :: b -> B.ByteString
    chromStart :: b -> Int
    chromEnd :: b -> Int

    toBed3 :: b -> BED3
    toBed3 bed = BED3 (chrom bed) (chromStart bed) (chromEnd bed)

    -- | return the size of a bed region
    size :: b -> Int
    size bed = chromEnd bed - chromStart bed

    hReadBed :: MonadIO m => Handle -> Source m b
    hReadBed h = do eof <- liftIO $ hIsEOF h
                    unless eof $ do
                        line <- liftIO $ B.hGetLine h
                        yield $ fromLine line
                        hReadBed h
    {-# INLINE hReadBed #-}

    -- | non-streaming version
    hReadBed' :: MonadIO m => Handle -> m [b]
    hReadBed' h = hReadBed h $$ sinkList
    {-# INLINE hReadBed' #-}

    readBed :: MonadIO m => FilePath -> Source m b
    readBed fl = do handle <- liftIO $ openFile fl ReadMode
                    hReadBed handle
                    liftIO $ hClose handle
    {-# INLINE readBed #-}

    -- | non-streaming version
    readBed' :: MonadIO m => FilePath -> m [b]
    readBed' fl = readBed fl $$ sinkList
    {-# INLINE readBed' #-}

    hWriteBed :: MonadIO m => Handle -> Sink b m ()
    hWriteBed handle = do
        x <- await
        case x of
            Nothing -> return ()
            Just bed -> (liftIO . B.hPutStrLn handle . toLine) bed >> hWriteBed handle
    {-# INLINE hWriteBed #-}

    hWriteBed' :: MonadIO m => Handle -> [b] -> m ()
    hWriteBed' handle beds = yieldMany beds $$ hWriteBed handle
    {-# INLINE hWriteBed' #-}

    writeBed :: MonadIO m => FilePath -> Sink b m ()
    writeBed fl = do handle <- liftIO $ openFile fl WriteMode
                     hWriteBed handle
                     liftIO $ hClose handle
    {-# INLINE writeBed #-}

    writeBed' :: MonadIO m => FilePath -> [b] -> m ()
    writeBed' fl beds = yieldMany beds $$ writeBed fl
    {-# INLINE writeBed' #-}

    {-# MINIMAL asBed, fromLine, toLine, chrom, chromStart, chromEnd #-}

type BEDTree a = M.HashMap B.ByteString (IM.IntervalMap Int a)

-- | convert a set of bed records to interval tree, with combining function for
-- equal keys
sortedBedToTree :: (BEDLike b, F.Foldable f)
                => (a -> a -> a)
                -> Sorted (f (b, a))
                -> BEDTree a
sortedBedToTree f (Sorted xs) =
      M.fromList
    . map ((head *** IM.fromAscListWith f) . unzip)
    . groupBy ((==) `on` fst)
    . map (\(bed, x) -> (chrom bed, (IM.IntervalCO (chromStart bed) (chromEnd bed), x)))
    . F.toList
    $ xs
{-# INLINE sortedBedToTree #-}

bedToTree :: BEDLike b
          => (a -> a -> a)
          -> [(b, a)]
          -> BEDTree a
bedToTree f xs =
      M.fromList
    . map ((head *** IM.fromAscListWith f) . unzip)
    . groupBy ((==) `on` fst)
    . map (\(bed, x) -> (chrom bed, (IM.IntervalCO (chromStart bed) (chromEnd bed), x)))
    . V.toList
    $ xs'
  where
    xs' = V.create $ do
        v <- V.unsafeThaw . V.fromList $ xs
        I.sortBy (compareBed `on` fst) v
        return v
{-# INLINE bedToTree #-}

-- | split a bed region into k consecutive subregions, discarding leftovers
splitBed :: BEDLike b => Int -> b -> [b]
splitBed k bed = map (uncurry (asBed chr)) . bins k $ (s, e)
  where
    chr = chrom bed
    s = chromStart bed
    e = chromEnd bed
{-# INLINE splitBed #-}

-- | split a bed region into consecutive fixed size subregions, discarding leftovers
splitBedBySize :: BEDLike b => Int -> b -> [b]
splitBedBySize k bed = map (uncurry (asBed chr)) . binBySize k $ (s, e)
  where
    chr = chrom bed
    s = chromStart bed
    e = chromEnd bed
{-# INLINE splitBedBySize #-}

-- | split a bed region into consecutive fixed size subregions, including leftovers
splitBedBySizeLeft :: BEDLike b => Int -> b -> [b]
splitBedBySizeLeft k bed = map (uncurry (asBed chr)) . binBySizeLeft k $ (s, e)
  where
    chr = chrom bed
    s = chromStart bed
    e = chromEnd bed
{-# INLINE splitBedBySizeLeft #-}

splitBedBySizeOverlap :: BEDLike b
                      => Int     -- ^ bin size
                      -> Int     -- ^ overlap size
                      -> b -> [b]
splitBedBySizeOverlap k o bed = map (uncurry (asBed chr)) .
    binBySizeOverlap k o $ (s, e)
  where
    chr = chrom bed
    s = chromStart bed
    e = chromEnd bed
{-# INLINE splitBedBySizeOverlap #-}

-- | a type to imply that underlying data structure is sorted
newtype Sorted b = Sorted {fromSorted :: b} deriving (Show, Read, Eq)

compareBed :: (BEDLike b1, BEDLike b2) => b1 -> b2 -> Ordering
compareBed x y = compare x' y'
  where
    x' = (chrom x, chromStart x, chromEnd x)
    y' = (chrom y, chromStart y, chromEnd y)
{-# INLINE compareBed #-}

-- | sort BED, first by chromosome (alphabetical order), then by chromStart, last by chromEnd
sortBed :: BEDLike b => [b] -> Sorted (V.Vector b)
sortBed beds = Sorted $ V.create $ do
    v <- V.unsafeThaw . V.fromList $ beds
    I.sortBy compareBed v
    return v
{-# INLINE sortBed #-}

-- | return records in A that are overlapped with records in B
intersectBed :: (BEDLike b1, BEDLike b2) => [b1] -> [b2] -> [b1]
intersectBed a b = intersectSortedBed a b'
  where
    b' = sortBed b
{-# INLINE intersectBed #-}

-- | return records in A that are overlapped with records in B
intersectSortedBed :: (BEDLike b1, BEDLike b2)
                   => [b1] -> Sorted (V.Vector b2) -> [b1]
intersectSortedBed a (Sorted b) = filter (not . IM.null . f) a
  where
    f bed = let chr = chrom bed
                interval = IM.IntervalCO (chromStart bed) $ chromEnd bed
            in IM.intersecting (M.lookupDefault IM.empty chr tree) interval
    tree = sortedBedToTree (\_ _ -> False) . Sorted . zip (V.toList b) . repeat $ False
{-# INLINE intersectSortedBed #-}

intersectBedWith :: (BEDLike b1, BEDLike b2)
                 => ([b2] -> a)
                 -> [b1]
                 -> [b2]
                 -> [(b1, a)]
intersectBedWith fn a = intersectSortedBedWith fn a . sortBed
{-# INLINE intersectBedWith #-}

intersectSortedBedWith :: (BEDLike b1, BEDLike b2)
                       => ([b2] -> a)
                       -> [b1]
                       -> Sorted (V.Vector b2)
                       -> [(b1, a)]
intersectSortedBedWith fn a (Sorted b) = map f a
  where
    f bed = let chr = chrom bed
                interval = IM.IntervalCO (chromStart bed) $ chromEnd bed
            in (bed, fn $ IM.elems $ IM.intersecting (M.lookupDefault IM.empty chr tree) interval)
    tree = sortedBedToTree const . Sorted . V.toList . V.zip b $ b
{-# INLINE intersectSortedBedWith #-}

isOverlapped :: (BEDLike b1, BEDLike b2) => b1 -> b2 -> Bool
isOverlapped bed1 bed2 = chr == chr' && not (e <= s' || e' <= s)
  where
    chr = chrom bed1
    s = chromStart bed1
    e = chromEnd bed1
    chr' = chrom bed2
    s' = chromStart bed2
    e' = chromEnd bed2

mergeBed :: (BEDLike b, Monad m) => [b] -> Source m b
mergeBed = mergeSortedBed . sortBed
{-# INLINE mergeBed #-}

mergeBedWith :: (BEDLike b, Monad m)
              => ([b] -> a) -> [b] -> Source m a
mergeBedWith f = mergeSortedBedWith f . sortBed
{-# INLINE mergeBedWith #-}

mergeSortedBed :: (BEDLike b, Monad m) => Sorted (V.Vector b) -> Source m b
mergeSortedBed = mergeSortedBedWith f
  where
    f xs = asBed (chrom $ head xs) lo hi
      where
        lo = minimum . map chromStart $ xs
        hi = maximum . map chromEnd $ xs
{-# INLINE mergeSortedBed #-}

mergeSortedBedWith :: (BEDLike b, Monad m)
                   => ([b] -> a) -> Sorted (V.Vector b) -> Source m a
mergeSortedBedWith mergeFn (Sorted beds) = do
    (_, r) <- V.foldM' f acc0 . V.tail $ beds
    yield $ mergeFn r
  where
    x0 = V.head beds
    acc0 = ((chrom x0, chromStart x0, chromEnd x0), [x0])
    f ((chr,lo,hi), acc) bed
        | chr /= chr' || s' > hi = yield (mergeFn acc) >>
                                   return ((chr',s',e'), [bed])
        | e' > hi = return ((chr',lo,e'), bed:acc)
        | otherwise = return ((chr,lo,hi), bed:acc)
      where
        chr' = chrom bed
        s' = chromStart bed
        e' = chromEnd bed
{-# INLINE mergeSortedBedWith #-}

-- * BED6 format

-- | BED6 format, as described in http://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
data BED = BED
    { _chrom :: !B.ByteString
    , _chromStart :: {-# UNPACK #-} !Int
    , _chromEnd :: {-# UNPACK #-} !Int
    , _name :: !(Maybe B.ByteString)
    , _score :: !(Maybe Double)
    , _strand :: !(Maybe Bool)  -- ^ True: "+", False: "-"
    } deriving (Eq, Show, Read)

instance Default BED where
    def = BED
        { _chrom = ""
        , _chromStart = 0
        , _chromEnd = 0
        , _name = Nothing
        , _score = Nothing
        , _strand = Nothing
        }

instance BEDLike BED where
    asBed chr s e = BED chr s e Nothing Nothing Nothing

    fromLine l = evalState (f (B.split '\t' l)) 1
      where
        f :: [B.ByteString] -> State Int BED
        f [] = do i <- get
                  if i <= 3 then error "Read BED fail: Incorrect number of fields"
                            else return def
        f (x:xs) = do
            i <- get
            put (i+1)
            bed <- f xs
            case i of
                1 -> return $ bed {_chrom = x}
                2 -> return $ bed {_chromStart = readInt x}
                3 -> return $ bed {_chromEnd = readInt x}
                4 -> return $ bed {_name = guard' x}
                5 -> return $ bed {_score = getScore x}
                6 -> return $ bed {_strand = getStrand x}
                _ -> return def

        guard' x | x == "." = Nothing
                 | otherwise = Just x
        getScore x | x == "." = Nothing
                   | otherwise = Just . readDouble $ x
        getStrand str | str == "-" = Just False
                      | str == "+" = Just True
                      | otherwise = Nothing
    {-# INLINE fromLine #-}

    toLine (BED f1 f2 f3 f4 f5 f6) = B.intercalate "\t" [ f1
                                                        , (B.pack.show) f2
                                                        , (B.pack.show) f3
                                                        , fromMaybe "." f4
                                                        , score'
                                                        , strand'
                                                        ]
      where
        strand' | f6 == Just True = "+"
                | f6 == Just False = "-"
                | otherwise = "."
        score' = case f5 of
                     Just x -> (B.pack.show) x
                     _ -> "."
    {-# INLINE toLine #-}

    chrom = _chrom
    chromStart = _chromStart
    chromEnd = _chromEnd

-- | retreive sequences
fetchSeq :: (BioSeq DNA a, MonadIO m) => Genome -> Conduit BED m (Either String (DNA a))
fetchSeq g = mapMC f
  where
    f (BED chr start end _ _ isForward) = do
        dna <- liftIO $ getSeq g (chr, start, end)
        return $ case isForward of
            Just False -> rc <$> dna
            _ -> dna
{-# INLINE fetchSeq #-}

fetchSeq' :: (BioSeq DNA a, MonadIO m) => Genome -> [BED] -> m [Either String (DNA a)]
fetchSeq' g beds = yieldMany beds $= fetchSeq g $$ sinkList
{-# INLINE fetchSeq' #-}

-- * BED3 format

data BED3 = BED3 !B.ByteString !Int !Int deriving (Eq, Show, Read)

instance Default BED3 where
    def = BED3 "" 0 0

instance BEDLike BED3 where
    asBed = BED3

    fromLine l = case B.split '\t' l of
                    (a:b:c:_) -> BED3 a (readInt b) $ readInt c
                    _ -> error "Read BED fail: Incorrect number of fields"
    {-# INLINE fromLine #-}

    toLine (BED3 f1 f2 f3) = B.intercalate "\t" [f1, (B.pack.show) f2, (B.pack.show) f3]
    {-# INLINE toLine #-}

    chrom (BED3 f1 _ _) = f1
    chromStart (BED3 _ f2 _) = f2
    chromEnd (BED3 _ _ f3) = f3

convert :: (BEDLike b1, BEDLike b2) => b1 -> b2
convert b = asBed (chrom b) (chromStart b) (chromEnd b)
{-# INLINE convert #-}
