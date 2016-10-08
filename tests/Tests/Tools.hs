module Tests.Tools (tests) where

import Bio.Utils.Functions
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Matrix.Unboxed as MU

tests :: TestTree
tests = testGroup "Test: Bio.Utils"
    [ testCase "quantileNormalization" quantileNormalizationTest
    ]

quantileNormalizationTest :: Assertion
quantileNormalizationTest = after @=? quantileNormalization before
  where
    before = MU.fromLists [ [2, 4, 4, 5]
                          , [5, 14, 4, 7]
                          , [4, 8, 6, 9]
                          , [3, 8, 5, 8]
                          , [3, 9, 3, 5] ]
    after = MU.fromLists [ [3.5, 3.5, 5.0, 5.0]
                         , [8.5, 8.5, 5.5, 5.5]
                         , [6.5, 5.0, 8.5, 8.5]
                         , [5.0, 5.5, 6.5, 6.5]
                         , [5.5, 6.5, 3.5, 3.5] ]
