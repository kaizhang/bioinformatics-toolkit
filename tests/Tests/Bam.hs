{-# LANGUAGE OverloadedStrings #-}

module Tests.Bam (tests) where

import           Bio.Data.Bam
import           Bio.Data.Bed
import           Conduit
import           Test.Tasty
import           Test.Tasty.HUnit

tests :: TestTree
tests = testGroup "Test: Bio.Data.Bam"
    [ testCase "bamToBed" bamToBedTest
    ]

bamToBedTest :: Assertion
bamToBedTest = do
    bed <- readBed' "tests/data/example.bed"
    bed' <- runBam $ readBam "tests/data/example.bam" =$= bamToBed $$ sinkList
    bed @=? bed'
