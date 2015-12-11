{-# LANGUAGE OverloadedStrings #-}

module Tests.Bed (tests) where

import Bio.Data.Bed
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Vector as V

tests :: TestTree
tests = testGroup "Test: Bio.Data.Bed"
    [ testCase "sortBed" sortBedTest
    , testCase "split" splitBedTest
    ]

sortBedTest :: Assertion
sortBedTest = do beds <- readBed' "tests/data/peaks.bed" :: IO [BED]
                 let (Sorted actual) = sortBed beds
                 expect <- fmap V.fromList $ readBed' "tests/data/peaks.sorted.bed"
                 expect @=? actual

splitBedTest :: Assertion
splitBedTest = (s1', s2', s3') @=? (s1, s2, s3)
  where
    bed = BED3 "chr1" 0 99
    s1 = splitBedBySize 20 bed
    s1' = map f [(0, 20), (20, 40), (40, 60), (60, 80)]
    s2 = splitBedBySizeLeft 20 bed
    s2' = map f [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]
    s3 = splitBedBySizeOverlap 20 10 bed
    s3' = map f [ ( 0, 20), (10, 30), (20, 40), (30, 50), (40, 60)
                , (50, 70), (60, 80), (70, 90), (80, 100) ]
    f (a,b) = asBed "chr1" a b
