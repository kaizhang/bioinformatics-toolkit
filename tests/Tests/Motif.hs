module Tests.Motif (tests) where

import Bio.Motif
import Test.Tasty
import Test.Tasty.HUnit
import Data.Packed.Matrix

tests :: TestTree
tests = testGroup "Test: Bio.Motif"
    [ testCase "IUPAC converting" toIUPACTest]

toIUPACTest :: Assertion
toIUPACTest = assertEqual "toIUPAC check" expect actual
  where
    expect = "SAA"
    actual = show . toIUPAC $ pwm
    pwm = PWM Nothing $ fromLists [ [0.1, 0.5, 0.4, 0]
                                  , [0.9, 0.1, 0, 0]
                                  , [0.51, 0.2, 0.09, 0.2]
                                  ]



