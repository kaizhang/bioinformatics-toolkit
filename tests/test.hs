{-# LANGUAGE OverloadedStrings #-}
import qualified Tests.Motif as Motif
import qualified Tests.Bed as Bed
import Test.Tasty
import Test.Tasty.Golden

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ Motif.tests
    , Bed.tests
    ]
