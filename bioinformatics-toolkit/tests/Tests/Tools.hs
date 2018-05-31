module Tests.Tools (tests) where

import Bio.Utils.Functions
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Matrix as MU
import Text.Printf (printf)

tests :: TestTree
tests = testGroup "Test: Bio.Utils"
    [ quantileNormalizationTest
    ]

quantileNormalizationTest :: TestTree
quantileNormalizationTest = testGroup "quantile normalization"
    [ testCase "case 1" $ y1 @=? quantileNormalization x1
    , testCase "case 2" $ MU.map (printf "%.2f") y2 @=?
        (MU.map (printf "%.2f") (quantileNormalization x2) :: MU.Matrix String)
    ]
  where
    x1 :: MU.Matrix Double
    x1 = MU.fromLists
        [ [2,  4, 4, 5]
        , [5, 14, 4, 7]
        , [4,  8, 6, 9]
        , [3,  8, 5, 8]
        , [3,  9, 3, 5] ]
    y1 :: MU.Matrix Double
    y1 = MU.fromLists
        [ [ 3.5,  3.5, 5.25, 4.25]
        , [ 8.5,  8.5, 5.25,  5.5]
        , [ 6.5, 5.25,  8.5,  8.5]
        , [5.25, 5.25,  6.5,  6.5]
        , [5.25,  6.5,  3.5, 4.25] ]
    x2 :: MU.Matrix Double
    x2 = MU.fromLists
        [ [5, 4, 3]
        , [2, 1, 4]
        , [3, 4, 6]
        , [4, 2, 8] ]
    y2 :: MU.Matrix Double
    y2 = MU.fromLists
        [ [ 5.666667,  5.166667, 2 ]
        , [ 2,  2, 3]
        , [3, 5.166667,  4.666667]
        , [4.666667,  3,  5.666667] ]
