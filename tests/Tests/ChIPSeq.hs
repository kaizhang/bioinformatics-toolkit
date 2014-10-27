{-# LANGUAGE OverloadedStrings #-}

module Tests.ChIPSeq (tests) where

import Bio.Data.Bed
import Bio.ChIPSeq
import Data.Conduit
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Vector as V
import System.IO.Unsafe (unsafePerformIO)

peaks :: IO [BED3]
peaks = readBed' "tests/data/peaks.bed"

tags :: Source IO BED
tags = readBed "tests/data/example.bed"

tests :: TestTree
tests = testGroup "Test: Bio.ChIPSeq"
    [ testCase "rpkm" testRPKM
    ]

testRPKM :: Assertion
testRPKM = a @=? b
  where
    (a, b) = unsafePerformIO $ do
        regions <- peaks
        r1 <- tags $$ rpkm regions
        r2 <- rpkmFromBam regions "tests/data/example.bam"
        return (V.toList r1, r2)
