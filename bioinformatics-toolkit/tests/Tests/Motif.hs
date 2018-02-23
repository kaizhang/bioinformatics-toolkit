{-# LANGUAGE OverloadedStrings #-}

module Tests.Motif (tests) where

import qualified Data.ByteString.Char8 as B
import           Data.Conduit
import qualified Data.Conduit.List     as CL
import           Data.Default.Class
import           System.Random
import           Test.Tasty
import           Test.Tasty.HUnit

import           Bio.Data.Fasta
import           Bio.Motif
import           Bio.Motif.Search
import           Bio.Seq

import qualified Data.Vector.Unboxed as U

dna :: DNA Basic
dna = case fromBS (B.pack $ map f $ take 5000 $ randomRs (0, 3) (mkStdGen 2)) of
    Left msg -> error msg
    Right x -> x
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
    , testCase "Max matching score" maxScTest
    , testCase "pValue calculation" pValueTest
    , testCase "CDF truncate test" cdfTruncateTest
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
    expect <- findTFBS def pwm dna (0.6 * optimalScore def pwm) True $$ CL.consume
    actual <- findTFBSSlow def pwm dna (0.6 * optimalScore def pwm) $$ CL.consume
    assertEqual "findTFBS" expect actual

maxScTest :: Assertion
maxScTest = do
    ms <- motifs
    let expect = map (\x -> maximum $ scores def (_pwm x) dna) ms
        actual = map (\x -> maxMatchSc def (_pwm x) dna) ms
    assertEqual "maxMatchSc" expect actual

pValueTest :: Assertion
pValueTest = do
    ms <- motifs
    let expect = map (approx . pValueToScoreExact 1e-4 def . _pwm) ms
        actual = map (approx . pValueToScore 1e-4 def . _pwm) ms
    assertEqual "pValueToScore" expect actual
  where
    approx x = round $ 10 * x

cdfTruncateTest :: Assertion
cdfTruncateTest = do
    ms <- motifs
    let expect = map (pValueToScore 1e-4 def . _pwm) ms
        actual = map (pValueToScore' 1e-4 def . _pwm) ms
    assertEqual "CDF truncate" expect actual
  where
    pValueToScore' p bg pwm = cdf' (truncateCDF 0.999 $ scoreCDF bg pwm) $ 1 - p
