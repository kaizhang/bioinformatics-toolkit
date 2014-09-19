{-# LANGUAGE OverloadedStrings #-}

module Tests.Motif (tests) where

import Bio.Seq
import Bio.Motif
import Bio.Motif.Search
import Test.Tasty
import Test.Tasty.HUnit
import System.Random
import qualified Data.Conduit.List as CL
import qualified Data.ByteString.Char8 as B
import Data.Conduit
import Control.Monad.Identity
import Data.Default.Generics

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

pwm :: PWM
pwm = readPWM "0.3 0.3 0.3 0.1\n0 0.5 0 0.5\n0.1 0.2 0.5 0.3\n0.1 0.1 0.1 0.7\n0 0 0 1\n0.25 0.25 0.25 0.25\n0.1 0.1 0.3 0.5\n0.25 0.25 0 0.5\n0.1 0.1 0.7 0.1\n0 0 0 1"

motif :: Motif
motif = Motif "test" pwm

tests :: TestTree
tests = testGroup "Test: Bio.Motif"
    [ testCase "IUPAC converting" toIUPACTest
    , testCase "Motif scanning" findTFBSTest
    ]

toIUPACTest :: Assertion
toIUPACTest = assertEqual "toIUPAC check" expect actual
  where
    expect = "SAA"
    actual = show . toIUPAC $ pwm

findTFBSTest :: Assertion
findTFBSTest = assertEqual "TFBS scan check" expect actual
  where
    expect = runIdentity $ findTFBS def motif dna (0.6 * optimalScore def pwm) $$ CL.consume
    actual = runIdentity $ findTFBS' def motif dna (0.6 * optimalScore def pwm) $$ CL.consume

{-
motifScoringCheck :: Assertion
motifScoringCheck = assertEqual "motif scoring check" expect actual
  where
    expect = 
    actual = 
-}
