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

    -- * BED6 format
    , BED(..)
    -- * BED3 format
    , BED3(..)
    -- * NarrowPeak format
    , NarrowPeak(..)

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

import Control.Arrow ((***))
import Control.Monad.State.Strict
import qualified Data.ByteString.Char8 as B
import Conduit
import Data.Default.Class (Default(..))
import Data.Function (on)
import qualified Data.Foldable as F
import qualified Data.HashMap.Strict as M
import qualified Data.IntervalMap.Strict as IM
import Data.List (groupBy, sortBy)
import Data.Maybe (fromMaybe, fromJust)
import qualified Data.Vector as V
import qualified Data.Vector.Algorithms.Intro as I
import System.IO
import Data.ByteString.Lex.Integral (packDecimal)
import Data.Double.Conversion.ByteString (toShortest)
import Data.Ord (comparing)

import Bio.Motif hiding (_name)
import Bio.Motif.Search
import Bio.Seq
import Bio.Seq.IO
import Bio.Utils.Misc ( binBySizeLeft, binBySize, binBySizeOverlap, bins
                      , readInt, readDouble )

-- | A class representing BED-like data, e.g., BED3, BED6 and BED12. BED format
-- uses 0-based index (see documentation).
class BEDLike b where
    -- | Construct bed record from chromsomoe, start location and end location
    asBed :: B.ByteString -> Int -> Int -> b

    -- | Convert bytestring to bed format
    fromLine :: B.ByteString -> b

    -- | Convert bed to bytestring
    toLine :: b -> B.ByteString

    -- | Field accessor
    chrom :: b -> B.ByteString
    chromStart :: b -> Int
    chromEnd :: b -> Int
    bedName :: b -> Maybe B.ByteString
    bedScore :: b -> Maybe Double
    bedStrand :: b -> Maybe Bool

    convert :: BEDLike b' => b' -> b
    convert bed = asBed (chrom bed) (chromStart bed) (chromEnd bed)
    {-# INLINE convert #-}

    -- | Return the size of a bed region.
    bedSize :: b -> Int
    bedSize bed = chromEnd bed - chromStart bed

    {-# MINIMAL asBed, fromLine, toLine, chrom, chromStart, chromEnd,
                bedName, bedScore, bedStrand #-}

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

    toLine (BED f1 f2 f3 f4 f5 f6) = B.intercalate "\t"
        [ f1, (B.pack.show) f2, (B.pack.show) f3, fromMaybe "." f4, score'
        , strand' ]
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
    bedName = _name
    bedScore = _score
    bedStrand = _strand

    convert bed = BED (chrom bed) (chromStart bed) (chromEnd bed) (bedName bed)
                      (bedScore bed) (bedStrand bed)

-- * BED3 format

data BED3 = BED3 !B.ByteString !Int !Int deriving (Eq, Show, Read)

instance BEDLike BED3 where
    asBed = BED3

    fromLine l = case B.split '\t' l of
                    (a:b:c:_) -> BED3 a (readInt b) $ readInt c
                    _ -> error "Read BED fail: Incorrect number of fields"
    {-# INLINE fromLine #-}

    toLine (BED3 a b c) = B.intercalate "\t"
        [a, fromJust $ packDecimal b, fromJust $ packDecimal c]
    {-# INLINE toLine #-}

    chrom (BED3 f1 _ _) = f1
    chromStart (BED3 _ f2 _) = f2
    chromEnd (BED3 _ _ f3) = f3
    bedName = const Nothing
    bedScore = const Nothing
    bedStrand = const Nothing

-- | ENCODE narrowPeak format: https://genome.ucsc.edu/FAQ/FAQformat.html#format12
data NarrowPeak = NarrowPeak
    { _npChrom :: !B.ByteString
    , _npStart :: !Int
    , _npEnd :: !Int
    , _npName :: !(Maybe B.ByteString)
    , _npScore :: !Double
    , _npStrand :: !(Maybe Bool)
    , _npSigal :: !Double
    , _npPvalue :: !(Maybe Double)
    , _npQvalue :: !(Maybe Double)
    , _npPeak :: !(Maybe Int)
    } deriving (Eq, Show, Read)

instance BEDLike NarrowPeak where
    asBed chr s e = NarrowPeak chr s e Nothing 0 Nothing 0 Nothing Nothing Nothing

    fromLine l = NarrowPeak a (readInt b) (readInt c)
        (if d == "." then Nothing else Just d)
        (readDouble e)
        (if f == "." then Nothing else if f == "+" then Just True else Just False)
        (readDouble g)
        (if readDouble h < 0 then Nothing else Just $ readDouble h)
        (if readDouble i < 0 then Nothing else Just $ readDouble i)
        (if readInt j < 0 then Nothing else Just $ readInt j)
      where
        (a:b:c:d:e:f:g:h:i:j:_) = B.split '\t' l
    {-# INLINE fromLine #-}

    toLine (NarrowPeak a b c d e f g h i j) = B.intercalate "\t"
        [ a, fromJust $ packDecimal b, fromJust $ packDecimal c, fromMaybe "." d
        , toShortest e
        , case f of
            Nothing -> "."
            Just True -> "+"
            _ -> "-"
        , toShortest g, fromMaybe "-1" $ fmap toShortest h
        , fromMaybe "-1" $ fmap toShortest i
        , fromMaybe "-1" $ fmap (fromJust . packDecimal) j
        ]
    {-# INLINE toLine #-}

    chrom = _npChrom
    chromStart = _npStart
    chromEnd = _npEnd
    bedName = _npName
    bedScore = Just . _npScore
    bedStrand = _npStrand

    convert bed = NarrowPeak (chrom bed) (chromStart bed) (chromEnd bed) (bedName bed)
        (fromMaybe 0 $ bedScore bed) (bedStrand bed) 0 Nothing Nothing Nothing

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

intersecting :: BEDLike b => BEDTree a -> b -> IM.IntervalMap Int a
intersecting tree x = IM.intersecting (M.lookupDefault IM.empty chr tree) interval
  where
    chr = chrom x
    interval = IM.IntervalCO (chromStart x) $ chromEnd x
{-# INLINE intersecting #-}

isIntersected :: BEDLike b => BEDTree a -> b -> Bool
isIntersected tree = not . IM.null . intersecting tree
{-# INLINE isIntersected #-}

sizeOverlapped :: (BEDLike b1, BEDLike b2) => b1 -> b2 -> Int
sizeOverlapped x y | chr1 /= chr2 = 0
                   | overlap < 0 = 0
                   | otherwise = overlap
  where
    overlap = minimum [e1 - s2, e2 - s1, e1 - s1, e2 - s2]
    chr1 = chrom x
    s1 = chromStart x
    e1 = chromEnd x
    chr2 = chrom y
    s2 = chromStart y
    e2 = chromEnd y


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
intersectBed :: (BEDLike b1, BEDLike b2, Monad m) => [b2] -> Conduit b1 m b1
intersectBed b = intersectSortedBed b'
  where
    b' = sortBed b
{-# INLINE intersectBed #-}

-- | return records in A that are overlapped with records in B
intersectSortedBed :: (BEDLike b1, BEDLike b2, Monad m)
                   => Sorted (V.Vector b2) -> Conduit b1 m b1
intersectSortedBed (Sorted b) = filterC (not . IM.null . intersecting tree)
  where
    tree = sortedBedToTree (\_ _ -> ()) . Sorted $ V.map (\x -> (x,())) b
{-# INLINE intersectSortedBed #-}

intersectBedWith :: (BEDLike b1, BEDLike b2, Monad m)
                 => ([b2] -> a)
                 -> [b2]
                 -> Conduit b1 m (b1, a)
intersectBedWith fn = intersectSortedBedWith fn . sortBed
{-# INLINE intersectBedWith #-}

intersectSortedBedWith :: (BEDLike b1, BEDLike b2, Monad m)
                       => ([b2] -> a)
                       -> Sorted (V.Vector b2)
                       -> Conduit b1 m (b1, a)
intersectSortedBedWith fn (Sorted b) = mapC f
  where
    f bed = (bed, fn $ IM.elems $ intersecting tree bed)
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
mergeSortedBedWith mergeFn (Sorted beds)
    | V.null beds = return ()
    | otherwise = do
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

-- | Split overlapped regions into non-overlapped regions. The input must be overlapped.
-- This function is usually used with `mergeBedWith`.
splitOverlapped :: BEDLike b => ([b] -> a) -> [b] -> [(BED3, a)]
splitOverlapped fun xs = filter ((>0) . bedSize . fst) $
    evalState (F.foldrM f [] $ init xs') x0
  where
    x0 = (\(a,b) -> (fromEither a, M.singleton (chromStart b, chromEnd b) b)) $ last xs'
    xs' = sortBy (comparing (fromEither . fst)) $ concatMap
        ( \x -> [(Left $ chromStart x, x), (Right $ chromEnd x, x)] ) xs
    f (i, x) acc = do
        (j, set) <- get
        let bed = (BED3 chr (fromEither i) j, fun $ M.elems set)
            set' = case i of
                Left _ -> M.delete (chromStart x, chromEnd x) set
                Right _ -> M.insert (chromStart x, chromEnd x) x set
        put (fromEither i, set')
        return (bed:acc)
    fromEither (Left x) = x
    fromEither (Right x) = x
    chr = chrom $ head xs
{-# INLINE splitOverlapped #-}

-- | Read records from a bed file handler in a streaming fashion.
hReadBed :: (BEDLike b, MonadIO m) => Handle -> Source m b
hReadBed h = do eof <- liftIO $ hIsEOF h
                unless eof $ do
                    line <- liftIO $ B.hGetLine h
                    yield $ fromLine line
                    hReadBed h
{-# INLINE hReadBed #-}

-- | Non-streaming version.
hReadBed' :: (BEDLike b, MonadIO m) => Handle -> m [b]
hReadBed' h = hReadBed h $$ sinkList
{-# INLINE hReadBed' #-}

-- | Read records from a bed file in a streaming fashion.
readBed :: (BEDLike b, MonadIO m) => FilePath -> Source m b
readBed fl = do handle <- liftIO $ openFile fl ReadMode
                hReadBed handle
                liftIO $ hClose handle
{-# INLINE readBed #-}

-- | Non-streaming version.
readBed' :: (BEDLike b, MonadIO m) => FilePath -> m [b]
readBed' fl = readBed fl $$ sinkList
{-# INLINE readBed' #-}

hWriteBed :: (BEDLike b, MonadIO m) => Handle -> Sink b m ()
hWriteBed handle = do
    x <- await
    case x of
        Nothing -> return ()
        Just bed -> (liftIO . B.hPutStrLn handle . toLine) bed >> hWriteBed handle
{-# INLINE hWriteBed #-}

hWriteBed' :: (BEDLike b, MonadIO m) => Handle -> [b] -> m ()
hWriteBed' handle beds = yieldMany beds $$ hWriteBed handle
{-# INLINE hWriteBed' #-}

writeBed :: (BEDLike b, MonadIO m) => FilePath -> Sink b m ()
writeBed fl = do handle <- liftIO $ openFile fl WriteMode
                 hWriteBed handle
                 liftIO $ hClose handle
{-# INLINE writeBed #-}

writeBed' :: (BEDLike b, MonadIO m) => FilePath -> [b] -> m ()
writeBed' fl beds = yieldMany beds $$ writeBed fl
{-# INLINE writeBed' #-}

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

-- | Identify motif binding sites
motifScan :: (BEDLike b, MonadIO m)
          => Genome -> [Motif] -> Bkgd -> Double -> Conduit b m BED
motifScan g motifs bg p = awaitForever $ \bed -> do
    let chr = chrom bed
        s = chromStart bed
        e = chromEnd bed
    r <- liftIO $ getSeq g (chr, s, e)
    case r of
        Left _ -> return ()
        Right dna -> mapM_ (getTFBS dna (chr, s)) motifs'
  where
    getTFBS dna (chr, s) (nm, (pwm, cutoff), (pwm', cutoff')) = toProducer
        ( (findTFBS bg pwm (dna :: DNA IUPAC) cutoff True =$=
            mapC (\i -> BED chr (s+i) (s+i+n) (Just nm) Nothing $ Just True)) >>
          (findTFBS bg pwm' dna cutoff' True =$=
            mapC (\i -> BED chr (s+i) (s+i+n) (Just nm) Nothing $ Just False)) )
      where
        n = Bio.Motif.size pwm
    motifs' = flip map motifs $ \(Motif nm pwm) ->
        let cutoff = pValueToScore p bg pwm
            cutoff' = pValueToScore p bg pwm'
            pwm' = rcPWM pwm
        in (nm, (pwm, cutoff), (pwm', cutoff'))
{-# INLINE motifScan #-}

-- | Retrieve motif matching scores
getMotifScore :: MonadIO m
              => Genome -> [Motif] -> Bkgd -> Conduit BED m BED
getMotifScore g motifs bg = awaitForever $ \(BED chr s e (Just nm) _ isForward) -> do
    r <- liftIO $ getSeq g (chr, s, e)
    let r' = case isForward of
            Just False -> rc <$> r
            _ -> r
    case r' of
        Left _ -> return ()
        Right dna -> do
            let pwm = M.lookupDefault (error "can't find motif with given name")
                      nm motifMap
                sc = score bg pwm (dna :: DNA IUPAC)
            yield $ BED chr s e (Just nm) (Just sc) isForward
  where
    motifMap = M.fromListWith (error "found motif with same name") $
        map (\(Motif nm pwm) -> (nm, pwm)) motifs
{-# INLINE getMotifScore #-}

getMotifPValue :: Monad m => [Motif] -> Bkgd -> Conduit BED m BED
getMotifPValue motifs bg = mapC $ \bed ->
    let nm = fromJust $ bedName bed
        sc = fromJust $ bedScore bed
        d = M.lookupDefault (error "can't find motif with given name")
                nm motifMap
        p = 1 - cdf d sc
     in bed{_score = Just p}
  where
    motifMap = M.fromListWith (error "getMotifPValue: found motif with same name") $
        map (\(Motif nm pwm) -> (nm, scoreCDF bg pwm)) motifs
{-# INLINE getMotifPValue #-}
