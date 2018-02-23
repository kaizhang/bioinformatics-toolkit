{-# LANGUAGE OverloadedStrings #-}
module Tests.Seq (tests) where

import           Bio.Seq
import           Data.Either
import qualified Data.HashMap.Strict as M
import           Test.Tasty
import           Test.Tasty.HUnit

tests :: TestTree
tests = testGroup "Test: Bio.Seq"
    [ testCase "nucleotideFreq" testNuclFreq
    ]

testNuclFreq :: Assertion
testNuclFreq = do
    let dna = fromRight undefined $ fromBS "ACTTCCCGGGD" :: DNA IUPAC
        test = M.lookupDefault undefined 'C' $ nucleotideFreq dna
        expected = 4
    test @=? expected
