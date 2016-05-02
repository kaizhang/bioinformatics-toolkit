{-# LANGUAGE OverloadedStrings #-}
module Tests.GREAT (tests) where

import           Bio.Data.Bed
import           Bio.GO.GREAT
import           Data.List
import           Data.Ord

import           Test.Tasty
import           Test.Tasty.HUnit

testData = (genes, results)
  where
    genes = [ (("chr1", 1000, True), 0)
            , (("chr1", 6000, True), 1)
            , (("chr2", 6000, True), 2)
            , (("chr2", 2000000, False), 3)
            ]
    results = [ (BED3 "chr1" 0 2000, 0)
              , (BED3 "chr1" 1000 1007000, 1)
              , (BED3 "chr2" 0 1003000, 2)
              , (BED3 "chr2" 1003000 3005000, 3)
              ]

tests :: TestTree
tests = testGroup "Test: Bio.GO.GREAT"
    [ testCase "getRegulatoryDomains" testAssoc
    ]

testAssoc :: Assertion
testAssoc = do
    let result = sortBy (comparing snd) $ getRegulatoryDomains
             (BasalPlusExtension 5000 1000 1000000) $ fst testData
    snd testData @=? result
