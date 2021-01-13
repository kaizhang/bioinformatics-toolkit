{-# LANGUAGE OverloadedStrings #-}

module Tests.Fastq (tests) where

import           Bio.Data.Fastq
import Data.Conduit.Zlib (ungzip, multiple, gzip)
import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy.Char8 as BL
import Lens.Micro
import           Conduit
import           Test.Tasty
import           Test.Tasty.Golden
import           Test.Tasty.HUnit
import qualified Data.HashMap.Strict as M
import Data.Ord
import Data.Maybe

tests :: TestTree
tests = testGroup "Test: Bio.Data.Fastq"
    [ testCase "FASTQ IO" fastqIO
    ]

fastqIO :: Assertion
fastqIO = do
    a <- fmap B.concat $ runResourceT $ runConduit $
        sourceFile "tests/data/test.fastq.gz" .| multiple ungzip .| sinkList
    b <- fmap B.concat $ runResourceT $ runConduit $
        sourceFile "tests/data/test_wrap.fastq.gz" .| multiple ungzip .| sinkList
    let fq1 = B.concat $ replicate 1000 a
        fq2 = B.concat $ replicate 1000 b
        r1 = B.concat $ runIdentity $ runConduit $ yield fq1 .|
            parseFastqC .| mapC fastqToByteString .| unlinesAsciiC .| sinkList
        r2 = B.concat $ runIdentity $ runConduit $ yield fq2 .|
            parseFastqC .| mapC fastqToByteString .| unlinesAsciiC .| sinkList
    [True, True] @=? [r1 == fq1, r2 == fq1]