{-# LANGUAGE OverloadedStrings #-}

module Tests.Fastq (tests) where

import           Bio.Data.Fastq
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
    [ fastqIO
    ]

fastqIO :: TestTree
fastqIO = testGroup "FASTQ Read/Write Test"
    [ goldenVsString "Standard" gold $ fmap (BL.fromStrict . B.concat) $ runResourceT $ runConduit $
          streamFastq gold .| mapC fastqToByteString .| unlinesAsciiC .| sinkList
    , goldenVsString "Wrap" gold $ fmap (BL.fromStrict . B.concat) $ runResourceT $ runConduit $
          streamFastq wrap .| mapC fastqToByteString .| unlinesAsciiC .| sinkList
    ]
  where
    gold = "tests/data/test.fastq" 
    wrap = "tests/data/test_wrap.fastq" 
    output = "out.fastq"
