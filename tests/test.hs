{-# LANGUAGE OverloadedStrings #-}
import qualified Tests.Motif as Motif
import qualified Tests.ChIPSeq as ChIPSeq
import Test.Tasty
import Test.Tasty.Golden

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ Motif.tests
    , ChIPSeq.tests
    ]
