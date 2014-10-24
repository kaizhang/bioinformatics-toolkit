{-# LANGUAGE OverloadedStrings #-}

module Tests.Bed (tests) where

import Bio.Data.Bed
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Vector as V

beds :: V.Vector BED3
beds = V.fromList [ BED3 "chr11" 1000 20000
                  , BED3 "chr2" 20 100
                  , BED3 "chr1" 1000 100000
                  , BED3 "chr1" 1 100
                  ]

sortedBeds :: V.Vector BED3
sortedBeds = V.fromList [ BED3 "chr1" 1 100
                        , BED3 "chr1" 1000 100000
                        , BED3 "chr2" 20 100
                        , BED3 "chr11" 1000 20000
                        ]

fromSorted :: Sorted a -> a
fromSorted (Sorted x) = x

tests :: TestTree
tests = testGroup "Test: Bio.Data.Bed"
    [ testCase "sortBED" $ assertEqual "sortBED" sortedBeds $ (fromSorted . sortBed) beds
    ]
