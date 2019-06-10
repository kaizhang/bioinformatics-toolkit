{-# LANGUAGE OverloadedStrings #-}

module Tests.Fastq (tests) where

import           Bio.Data.Fastq
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
    [ fastqIO
    ]

fastqIO :: TestTree
fastqIO = goldenVsFile "FASTQ Read/Write Test" input output io
  where
    io = runResourceT $ runConduit $ streamFastq input .| sinkFastq output
    input = "tests/data/test.fastq" 
    output = "out.fastq"
