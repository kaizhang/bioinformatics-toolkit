{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

module Bio.Data.Bed.Types
    ( BEDLike(..)
    , BEDConvert(..)
    , BED(..)
    , BED3(..)
    , NarrowPeak(..)
    , npSignal
    , npPvalue
    , npQvalue
    , npPeak
    , BroadPeak(..)
    , bpSignal
    , bpPvalue
    , bpQvalue
    , BEDExt(..)
    , _bed
    , _data
    , BEDTree
    , Sorted(..)
    ) where

import Lens.Micro
import Lens.Micro.TH (makeLensesFor)
import qualified Data.ByteString.Char8             as B
import           Data.ByteString.Lex.Integral      (packDecimal)
import           Data.Default.Class                (Default (..))
import           Data.Double.Conversion.ByteString (toShortest)
import qualified Data.HashMap.Strict               as M
import qualified Data.IntervalMap.Strict           as IM
import           Data.Maybe                        (fromJust, fromMaybe)

import           Bio.Utils.Misc                    (readDouble, readInt)

readDoubleNonnegative :: B.ByteString -> Maybe Double
readDoubleNonnegative x | v < 0 = Nothing
                        | otherwise = Just v
  where
    v = readDouble x
{-# INLINE readDoubleNonnegative #-}

readIntNonnegative :: B.ByteString -> Maybe Int
readIntNonnegative x | v < 0 = Nothing
                     | otherwise = Just v
  where
    v = readInt x
{-# INLINE readIntNonnegative #-}

-- | A class representing BED-like data, e.g., BED3, BED6 and BED12. BED format
-- uses 0-based index (see documentation).
class BEDLike b where
    -- | Field lens
    chrom :: Lens' b B.ByteString
    chromStart :: Lens' b Int
    chromEnd :: Lens' b Int
    name :: Lens' b (Maybe B.ByteString)
    score :: Lens' b (Maybe Int)
    strand :: Lens' b (Maybe Bool)

    -- | Return the size of a bed region.
    size :: b -> Int
    size bed = bed^.chromEnd - bed^.chromStart
    {-# INLINE size #-}

    {-# MINIMAL chrom, chromStart, chromEnd, name, score, strand #-}

class BEDLike b => BEDConvert b where
    -- | Construct bed record from chromsomoe, start location and end location
    asBed :: B.ByteString -> Int -> Int -> b

    -- | Convert bytestring to bed format
    fromLine :: B.ByteString -> b

    -- | Convert bed to bytestring
    toLine :: b -> B.ByteString

    convert :: BEDLike b' => b' -> b
    convert bed = asBed (bed^.chrom) (bed^.chromStart) (bed^.chromEnd)
    {-# INLINE convert #-}

    {-# MINIMAL asBed, fromLine, toLine #-}

-- * BED6 format

-- | BED6 format, as described in http://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
data BED = BED
    { _bed_chrom      :: !B.ByteString
    , _bed_chromStart :: !Int
    , _bed_chromEnd   :: !Int
    , _bed_name       :: !(Maybe B.ByteString)
    , _bed_score      :: !(Maybe Int)
    , _bed_strand     :: !(Maybe Bool)  -- ^ True: "+", False: "-"
    } deriving (Eq, Show, Read)

instance Ord BED where
    compare (BED x1 x2 x3 x4 x5 x6) (BED y1 y2 y3 y4 y5 y6) =
        compare (x1,x2,x3,x4,x5,x6) (y1,y2,y3,y4,y5,y6)

instance BEDLike BED where
    chrom = lens _bed_chrom (\bed x -> bed { _bed_chrom = x })
    chromStart = lens _bed_chromStart (\bed x -> bed { _bed_chromStart = x })
    chromEnd = lens _bed_chromEnd (\bed x -> bed { _bed_chromEnd = x })
    name = lens _bed_name (\bed x -> bed { _bed_name = x })
    score = lens _bed_score (\bed x -> bed { _bed_score = x })
    strand = lens _bed_strand (\bed x -> bed { _bed_strand = x })

instance BEDConvert BED where
    asBed chr s e = BED chr s e Nothing Nothing Nothing

    fromLine l = f $ take 6 $ B.split '\t' l
      where
        f [f1,f2,f3,f4,f5,f6] = BED f1 (readInt f2) (readInt f3) (getName f4)
            (getScore f5) (getStrand f6)
        f [f1,f2,f3,f4,f5] = BED f1 (readInt f2) (readInt f3) (getName f4)
            (getScore f5) Nothing
        f [f1,f2,f3,f4] = BED f1 (readInt f2) (readInt f3) (getName f4)
            Nothing Nothing
        f [f1,f2,f3] = asBed f1 (readInt f2) (readInt f3)
        f _ = error "Read BED fail: Not enough fields!"
        getName x | x == "." = Nothing
                  | otherwise = Just x
        getScore x | x == "." = Nothing
                   | otherwise = readInt x
        getStrand str | str == "-" = Just False
                      | str == "+" = Just True
                      | otherwise = Nothing
    {-# INLINE fromLine #-}

    toLine (BED f1 f2 f3 f4 f5 f6) = B.intercalate "\t"
        [ f1, fromJust $ packDecimal f2, fromJust $ packDecimal f3
        , fromMaybe "." f4, score', strand' ]
      where
        strand' | f6 == Just True = "+"
                | f6 == Just False = "-"
                | otherwise = "."
        score' = case f5 of
                     Just x -> fromJust $ packDecimal x
                     _      -> "."
    {-# INLINE toLine #-}

    convert bed = BED (bed^.chrom) (bed^.chromStart) (bed^.chromEnd) (bed^.name)
                      (bed^.score) (bed^.strand)

-- * BED3 format

data BED3 = BED3
    { _bed3_chrom       :: !B.ByteString
    , _bed3_chrom_start :: !Int
    , _bed3_chrom_end   :: !Int
    } deriving (Eq, Show, Read)

instance Ord BED3 where
    compare (BED3 x1 x2 x3) (BED3 y1 y2 y3) = compare (x1,x2,x3) (y1,y2,y3)

instance BEDLike BED3 where
    chrom = lens _bed3_chrom (\bed x -> bed { _bed3_chrom = x })
    chromStart = lens _bed3_chrom_start (\bed x -> bed { _bed3_chrom_start = x })
    chromEnd = lens _bed3_chrom_end (\bed x -> bed { _bed3_chrom_end = x })
    name = lens (const Nothing) (\bed _ -> bed)
    score = lens (const Nothing) (\bed _ -> bed)
    strand = lens (const Nothing) (\bed _ -> bed)

instance BEDConvert BED3 where
    asBed = BED3

    fromLine l = case B.split '\t' l of
                    (a:b:c:_) -> BED3 a (readInt b) $ readInt c
                    _ -> error "Read BED fail: Incorrect number of fields"
    {-# INLINE fromLine #-}

    toLine (BED3 a b c) = B.intercalate "\t"
        [a, fromJust $ packDecimal b, fromJust $ packDecimal c]
    {-# INLINE toLine #-}

-- | ENCODE narrowPeak format: https://genome.ucsc.edu/FAQ/FAQformat.html#format12
data NarrowPeak = NarrowPeak
    { _npChrom  :: !B.ByteString
    , _npStart  :: !Int
    , _npEnd    :: !Int
    , _npName   :: !(Maybe B.ByteString)
    , _npScore  :: !Int
    , _npStrand :: !(Maybe Bool)
    , _npSignal  :: !Double
    , _npPvalue :: !(Maybe Double)
    , _npQvalue :: !(Maybe Double)
    , _npPeak   :: !(Maybe Int)
    } deriving (Eq, Show, Read)

makeLensesFor [ ("_npSignal", "npSignal")
              , ("_npPvalue", "npPvalue")
              , ("_npQvalue", "npQvalue")
              , ("_npPeak", "npPeak")
              ] ''NarrowPeak

instance BEDLike NarrowPeak where
    chrom = lens _npChrom (\bed x -> bed { _npChrom = x })
    chromStart = lens _npStart (\bed x -> bed { _npStart = x })
    chromEnd = lens _npEnd (\bed x -> bed { _npEnd = x })
    name = lens _npName (\bed x -> bed { _npName = x })
    score = lens (Just . _npScore) (\bed x -> bed { _npScore = fromJust x })
    strand = lens _npStrand (\bed x -> bed { _npStrand = x })

instance BEDConvert NarrowPeak where
    asBed chr s e = NarrowPeak chr s e Nothing 0 Nothing 0 Nothing Nothing Nothing

    fromLine = go . B.split '\t'
      where
        go [a,b,c] = convert $ BED3 a (readInt b) $ readInt c
        go (a:b:c:d:e:f:g:h:i:j:_) = NarrowPeak a (readInt b) (readInt c)
            (if d == "." then Nothing else Just d)
            (readInt e)
            (if f == "." then Nothing else if f == "+" then Just True else Just False)
            (readDouble g)
            (readDoubleNonnegative h)
            (readDoubleNonnegative i)
            (readIntNonnegative j)
        go x = error $ "Cannot parse line: " <> show x
    {-# INLINE fromLine #-}

    toLine (NarrowPeak a b c d e f g h i j) = B.intercalate "\t"
        [ a, fromJust $ packDecimal b, fromJust $ packDecimal c, fromMaybe "." d
        , fromJust $ packDecimal e
        , case f of
            Nothing   -> "."
            Just True -> "+"
            _         -> "-"
        , toShortest g, fromMaybe "-1" $ fmap toShortest h
        , fromMaybe "-1" $ fmap toShortest i
        , fromMaybe "-1" $ fmap (fromJust . packDecimal) j
        ]
    {-# INLINE toLine #-}

    convert bed = NarrowPeak (bed^.chrom) (bed^.chromStart) (bed^.chromEnd) (bed^.name)
        (fromMaybe 0 $ bed^.score) (bed^.strand) 0 Nothing Nothing Nothing

-- | ENCODE broadPeak format: https://genome.ucsc.edu/FAQ/FAQformat.html#format13
data BroadPeak = BroadPeak
    { _bpChrom  :: !B.ByteString
    , _bpStart  :: !Int
    , _bpEnd    :: !Int
    , _bpName   :: !(Maybe B.ByteString)
    , _bpScore  :: !Int
    , _bpStrand :: !(Maybe Bool)
    , _bpSignal  :: !Double
    , _bpPvalue :: !(Maybe Double)
    , _bpQvalue :: !(Maybe Double)
    } deriving (Eq, Show, Read)

makeLensesFor [ ("_bpSignal", "bpSignal")
              , ("_bpPvalue", "bpPvalue")
              , ("_bpQvalue", "bpQvalue")
              ] ''BroadPeak

instance BEDLike BroadPeak where
    chrom = lens _bpChrom (\bed x -> bed { _bpChrom = x })
    chromStart = lens _bpStart (\bed x -> bed { _bpStart = x })
    chromEnd = lens _bpEnd (\bed x -> bed { _bpEnd = x })
    name = lens _bpName (\bed x -> bed { _bpName = x })
    score = lens (Just . _bpScore) (\bed x -> bed { _bpScore = fromJust x })
    strand = lens _bpStrand (\bed x -> bed { _bpStrand = x })

instance BEDConvert BroadPeak where
    asBed chr s e = BroadPeak chr s e Nothing 0 Nothing 0 Nothing Nothing

    fromLine l = BroadPeak a (readInt b) (readInt c)
        (if d == "." then Nothing else Just d)
        (readInt e)
        (if f == "." then Nothing else if f == "+" then Just True else Just False)
        (readDouble g)
        (readDoubleNonnegative h)
        (readDoubleNonnegative i)
      where
        (a:b:c:d:e:f:g:h:i:_) = B.split '\t' l
    {-# INLINE fromLine #-}

    toLine (BroadPeak a b c d e f g h i) = B.intercalate "\t"
        [ a, fromJust $ packDecimal b, fromJust $ packDecimal c, fromMaybe "." d
        , fromJust $ packDecimal e
        , case f of
            Nothing   -> "."
            Just True -> "+"
            _         -> "-"
        , toShortest g, fromMaybe "-1" $ fmap toShortest h
        , fromMaybe "-1" $ fmap toShortest i
        ]
    {-# INLINE toLine #-}

    convert bed = BroadPeak (bed^.chrom) (bed^.chromStart) (bed^.chromEnd) (bed^.name)
        (fromMaybe 0 $ bed^.score) (bed^.strand) 0 Nothing Nothing


data BEDExt bed a = BEDExt
    { _ext_bed :: bed
    , _ext_data :: a
    } deriving (Eq, Show, Read)

makeLensesFor [("_ext_bed", "_bed"), ("_ext_data", "_data")] ''BEDExt

instance BEDLike bed => BEDLike (BEDExt bed a) where
    chrom = _bed . chrom
    chromStart = _bed . chromStart
    chromEnd = _bed . chromEnd
    name = _bed . name
    score = _bed . score
    strand = _bed . strand

instance (Default a, Read a, Show a, BEDConvert bed) => BEDConvert (BEDExt bed a) where
    asBed chr s e = BEDExt (asBed chr s e) def

    fromLine l = let (a, b) = B.breakEnd (=='\t') l
                 in BEDExt (fromLine $ B.init a) $ read $ B.unpack b
    {-# INLINE fromLine #-}

    toLine (BEDExt bed a) = toLine bed <> "\t" <> B.pack (show a)
    {-# INLINE toLine #-}

type BEDTree a = M.HashMap B.ByteString (IM.IntervalMap Int a)

-- | a type to imply that underlying data structure is sorted
newtype Sorted b = Sorted {fromSorted :: b} deriving (Show, Read, Eq)