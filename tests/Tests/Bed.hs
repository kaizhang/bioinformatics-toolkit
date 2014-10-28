{-# LANGUAGE OverloadedStrings #-}

module Tests.Bed (tests) where

import Bio.Data.Bed
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Vector as V

tests :: TestTree
tests = testGroup "Test: Bio.Data.Bed"
    [ testCase "sortBed" sortBedTest
    ]

sortBedTest :: Assertion
sortBedTest = do beds <- readBed' "tests/data/peaks.bed" :: IO [BED]
                 let (Sorted actual) = sortBed . V.fromList $ beds
                 expect <- fmap V.fromList $ readBed' "tests/data/peaks.sorted.bed"
                 expect @=? actual
