{-# LANGUAGE OverloadedStrings #-}
module Tests.GREAT (tests) where

import           Bio.Data.Bed
import           Bio.GO.GREAT
import           Conduit
import           Control.Monad.Identity (runIdentity)
import           Data.Function          (on)
import           Data.List
import           Data.Ord

import           Test.Tasty
import           Test.Tasty.HUnit

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

loops = [ (BED3 "chr1" 900 1200, BED3 "chr1" 10000 11000)
        , (BED3 "chr2" 0 100, BED3 "chr2" 5500 6500)
        ]

results3D = [ (BED3 "chr1" 0 2000, 0)
            , (BED3 "chr1" 1000 7000, 1)
            , (BED3 "chr2" 1000 7000, 2)
            , (BED3 "chr2" 1999000 2005000, 3)
            , (BED3 "chr1" 10000 11000, 0)
            , (BED3 "chr1" 10000 11000, 1)
            , (BED3 "chr2" 0 100, 2)
            ]

tests :: TestTree
tests = testGroup "Test: Bio.GO.GREAT"
    [ testCase "getRegulatoryDomains" testAssoc
    , testCase "get3DRegulatoryDomains" testAssoc3D
    ]

testAssoc :: Assertion
testAssoc = results @=? sortBy (comparing snd)
    (getRegulatoryDomains (BasalPlusExtension 5000 1000 1000000) genes)

testAssoc3D :: Assertion
testAssoc3D = sortBy (compareBed `on` fst) results3D @=?
    sortBy (compareBed `on` fst) (runIdentity $ yieldMany loops =$=
        get3DRegulatoryDomains genes 5000 1000 $$ sinkList)
