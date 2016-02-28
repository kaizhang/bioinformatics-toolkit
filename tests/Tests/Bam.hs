{-# LANGUAGE OverloadedStrings #-}

module Tests.Bam (tests) where

import Bio.Data.Bed
import Bio.Data.Bam
import Conduit
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Vector as V

tests :: TestTree
tests = testGroup "Test: Bio.Data.Bam"
    [ testCase "bamToBed" bamToBedTest
    ]

bamToBedTest :: Assertion
bamToBedTest = do
    bed <- readBed' "tests/data/example.bed"
    bed' <- readBam "tests/data/example.bam" =$= bamToBed $$ sinkList
    bed @=? bed'
