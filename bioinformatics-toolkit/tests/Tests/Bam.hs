{-# LANGUAGE OverloadedStrings #-}

module Tests.Bam (tests) where

import           Bio.Data.Bam
import           Bio.Data.Bed
import           Bio.Utils.Misc        (readInt)
import           Conduit
import           Control.Monad         (forM_)
import qualified Data.ByteString.Char8 as B
import           Data.Tuple            (swap)
import           Test.Tasty
import           Test.Tasty.Golden
import           Test.Tasty.HUnit

tests :: TestTree
tests = testGroup "Test: Bio.Data.Bam"
    [ bamIOTest
    , testCase "bamToBed" bamToBedTest
    , testCase "sortedBamToBedPE" sortedBamToBedPETest
    ]

bamIOTest :: TestTree
bamIOTest = do
    goldenVsFile "BAM Read/Write Test" input output io
  where
    io = do
        header <- getBamHeader input
        runResourceT $ runConduit $ streamBam input .| sinkBam output header
    input = "tests/data/example.bam"
    output = "tests/data/example_copy.bam"

bamToBedTest :: Assertion
bamToBedTest = do
    bed <- readBed' "tests/data/example.bed"
    bed' <- do
        let input = "tests/data/example.bam" 
        header <- getBamHeader input
        runResourceT $ runConduit $ streamBam input .| bamToBedC header .| sinkList
    forM_ (zip bed bed') $ \(a,b) -> 
        if a == b then return () else error $ show (a,b)
    (bed == bed') @? show (head bed, head bed')

sortedBamToBedPETest :: Assertion
sortedBamToBedPETest = do
    bedpe <- readBedPE "tests/data/pairedend.bedpe" :: IO [(BED3, BED3)]
    bedpe' <- do
        let input = "tests/data/pairedend.bam"
        header <- getBamHeader input
        runResourceT $ runConduit $ streamBam input .|
            sortedBamToBedPE header .|
            mapC (\(x,y) -> (convert x, convert y)) .| sinkList
    forM_ (zip bedpe bedpe') $ \(b1, b2) -> (b1 == b2 || b1 == swap b2) @? show (b1,b2)
  where
    readBedPE fl = do
        c <- B.readFile fl
        return $ map (f . B.split '\t') $ B.lines c
    f (f1:f2:f3:f4:f5:f6:_) = ( asBed f1 (readInt f2) (readInt f3)
                              , asBed f4 (readInt f5) (readInt f6) )
