{-# LANGUAGE OverloadedStrings #-}

module Tests.Motif (tests) where

import Bio.Seq
import Bio.Motif
import Bio.Motif.Search
import Bio.Data.Fasta
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.Golden
import System.Random
import qualified Data.Conduit.List as CL
import qualified Data.ByteString.Char8 as B
import Data.Conduit

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
motifs = readFasta "tests/data/motifs.fasta" $$ CL.consume

f :: IO ()
f = do m <- motifs
       writeFasta "out.fasta" m

tests :: TestTree
tests = testGroup "Test: Bio.Motif"
    [ --testCase "IUPAC converting" toIUPACTest
    --, testCase "Motif scanning" findTFBSTest
     goldenVsFile "Motif IO" "tests/data/motifs.fasta" "out.fasta" f
    ]

{-
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
    -}
