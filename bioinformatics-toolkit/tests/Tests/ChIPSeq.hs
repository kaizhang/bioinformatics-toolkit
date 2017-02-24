module Tests.ChIPSeq (tests) where

import Bio.Data.Bed
import Bio.ChIPSeq
import Data.Conduit
import qualified Data.Conduit.List as CL
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Vector as V
import Text.Printf

peaks :: IO [BED3]
peaks = readBed' "tests/data/peaks.bed"

tags :: Source IO BED
tags = readBed "tests/data/example.bed"

tests :: TestTree
tests = testGroup "Test: Bio.ChIPSeq"
    [ testCase "rpkm" testRPKM
    ]

testRPKM :: Assertion
testRPKM = do regions <- peaks
              r1 <- tags $$ rpkmBed regions
              r2 <- CL.sourceList regions $= rpkmBam "tests/data/example.bam" $$ CL.consume
              let r1' = map (printf "%0.6f") . V.toList $ r1 :: [String]
                  r2' = map (printf "%0.6f") r2
              r1' @=? r2'
