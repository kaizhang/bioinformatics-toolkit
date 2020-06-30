{-# LANGUAGE OverloadedStrings #-}

module Tests.Bed (tests) where

import           Bio.Data.Bed
import           Bio.Data.Bed.Types
import           Bio.Data.Bed.Utils
import Bio.Utils.BitVector hiding (size)
import Lens.Micro
import           Conduit
import           Data.Function    (on)
import           Data.List        (sortBy, sort)
import qualified Data.Vector      as V
import           Test.Tasty
import           Test.Tasty.HUnit
import qualified Data.HashMap.Strict as M
import Data.Ord
import Data.Maybe
import Bio.RealWorld.GENCODE
import Data.Conduit.Zlib 

tests :: TestTree
tests = testGroup "Test: Bio.Data.Bed"
    [ testCase "sortBed" sortBedTest
    , testCase "split" splitBedTest
    , testCase "splitOverlapped" splitOverlappedTest
    , testCase "mergeBed" mergeBedTest
    , testCase "intersectBed" intersectBedTest
    , testCase "baseMap" baseMapTest
    , testCase "domain assignment" geneDomainTest
    ]

sortBedTest :: Assertion
sortBedTest = do
    beds <- runResourceT $ runConduit $
        streamBedGzip "tests/data/peaks.bed.gz" .| sinkList :: IO [BED3]
    let (Sorted actual) = sortBed beds
    expect <- runResourceT $ runConduit $
        streamBedGzip "tests/data/peaks.sorted.bed.gz" .| sinkVector
    expect @=? actual

splitBedTest :: Assertion
splitBedTest = (s1', s2', s3') @=? (s1, s2, s3)
  where
    bed = asBed "chr1" 0 99
    s1 = splitBedBySize 20 bed
    s1' = map f [(0, 20), (20, 40), (40, 60), (60, 80)]
    s2 = splitBedBySizeLeft 20 bed
    s2' = map f [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]
    s3 = splitBedBySizeOverlap 20 10 bed
    s3' = map f [ ( 0, 20), (10, 30), (20, 40), (30, 50), (40, 60)
                , (50, 70), (60, 80), (70, 90), (80, 100) ]
    f (a,b) = asBed "chr1" a b :: BED3

splitOverlappedTest :: Assertion
splitOverlappedTest = expect @=? result
  where
    input :: [BED3]
    input = map (\(a,b) -> asBed "chr1" a b)
        [ (0, 100)
        , (10, 20)
        , (50, 150)
        , (120, 160)
        , (155, 200)
        , (155, 220)
        , (0, 10)
        , (0, 10)
        ]
    expect = map (\((a,b), x) -> (asBed "chr1" a b, x))
        [ ((0, 10), 3)
        , ((10, 20), 2)
        , ((20, 50), 1)
        , ((50, 100), 2)
        , ((100, 120), 1)
        , ((120, 150), 2)
        , ((150, 155), 1)
        , ((155, 160), 3)
        , ((160, 200), 2)
        , ((200, 220), 1)
        ]
    --result = sortBy (compareBed `on` fst) $ splitOverlapped length input
    result = sortBy (compareBed `on` fst) $ countOverlapped input

mergeBedTest :: Assertion
mergeBedTest = expect @=? result
  where
    input :: [BED3]
    input = map (\(a,b) -> asBed "chr1" a b)
        [ (0, 100)
        , (10, 20)
        , (50, 150)
        , (120, 160)
        , (155, 200)
        , (155, 220)
        , (500, 1000)
        ]
    expect = map (\(a,b) -> asBed "chr1" a b)
        [ (0, 220)
        , (500, 1000)
        ]
    result = runIdentity $ runConduit $ mergeBed input .| sinkList

intersectBedTest :: Assertion
intersectBedTest = do
    expect <- runResourceT $ runConduit $
        streamBedGzip "tests/data/example_intersect_peaks.bed.gz" .| sinkList :: IO [BED3]
    peaks <- runResourceT $ runConduit $
        streamBedGzip "tests/data/peaks.bed.gz" .| sinkList :: IO [BED3]
    result <- runResourceT $ runConduit $
        streamBedGzip "tests/data/example.bed.gz" .| intersectBed peaks .| sinkList
    expect @=? result

baseMapTest :: Assertion
baseMapTest = do
    BaseMap bv <- runResourceT $ runConduit $ streamBedGzip "tests/data/example.bed.gz" .|
        baseMap [("chr1", 300000000)]
    let res = M.lookupDefault undefined "chr1" $
            fmap (map fst . filter snd . zip [0..] . toList) bv
    expect <- runResourceT $ runConduit $ streamBedGzip "tests/data/example.bed.gz" .|
        concatMapC f .| sinkList
    sort expect @=? sort res
  where
    f :: BED -> Maybe Int
    f bed = if bed^.chrom == "chr1"
        then Just $ if fromJust (bed^.strand)
            then bed^.chromStart
            else bed^.chromEnd - 1
        else Nothing

geneDomainTest :: Assertion
geneDomainTest = do
    genes <- runResourceT $ runConduit $ sourceFile "tests/data/genes.gtf.gz" .|
        multiple ungzip .| readGenesC
    let promoters = bedToTree const $ zip (concatMap (getPromoters 5000 1000) genes) $ repeat ()
        domain = getDomains 1000000 $ concatMap (getPromoters 5000 1000) genes
        f x = not (isIntersected promoters x) && (c1 || c2)
          where
            c1 = isIntersected promoters (chromStart %~ (`subtract` 1) $ x) &&
                isIntersected promoters (chromEnd %~ (+1) $ x) &&
                size x < 1000000
            c2 = isIntersected promoters (chromStart %~ (`subtract` 1) $ chromEnd %~ (+1) $ x) &&
                size x == 1000000
    all f domain @=? True