{-# LANGUAGE OverloadedStrings #-}
import qualified Tests.Motif as Motif
import Test.Tasty
import Test.Tasty.Golden

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ Motif.tests ]
