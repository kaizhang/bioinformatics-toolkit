{-# LANGUAGE OverloadedStrings #-}

module Tests.Motif (tests) where

import Test.Tasty
import Test.Tasty.HUnit
import System.Random
import Data.Default.Class
import qualified Data.Conduit.List as CL
import qualified Data.ByteString.Char8 as B
import Data.Conduit

import Bio.Seq
import Bio.Motif
import Bio.Motif.Search
import Bio.Data.Fasta

dna :: DNA Basic
dna = fromBS $ B.pack $ map f $ take 5000 $ randomRs (0, 3) (mkStdGen 2)
  where
    f :: Int -> Char
    f x = case x of
        0 -> 'A'
        1 -> 'C'
        2 -> 'G'
        3 -> 'T'
        _ -> undefined

motifs :: IO [Motif]
motifs = readFasta' "tests/data/motifs.fasta"

tests :: TestTree
tests = testGroup "Test: Bio.Motif"
    [ --testCase "IUPAC converting" toIUPACTest
      testCase "TFBS scanning" findTFBSTest
    , testCase "pValue calculation" pValueTest
    ]


{-
toIUPACTest :: Assertion
toIUPACTest = assertEqual "toIUPAC check" expect actual
  where
    expect = "SAA"
    actual = show . toIUPAC $ pwm
    -}

findTFBSTest :: Assertion
findTFBSTest = do
    ms <- motifs
    let (Motif _ pwm) = head ms
    expect <- findTFBS def pwm dna (0.6 * optimalScore def pwm) $$ CL.consume
    actual <- findTFBSSlow def pwm dna (0.6 * optimalScore def pwm) $$ CL.consume
    assertEqual "findTFBS" expect actual

pValueTest :: Assertion
pValueTest = do
    ms <- motifs
    let expect = map (pValueToScoreExact 1e-4 def . _pwm) ms
        actual = map (pValueToScore 1e-4 def . _pwm) ms
    assertEqual "pValueToScore" expect actual
