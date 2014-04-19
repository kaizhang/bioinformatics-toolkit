{-# LANGUAGE OverloadedStrings, UnicodeSyntax, TemplateHaskell, BangPatterns #-}
import qualified Data.ByteString.Lazy.Char8 as B
import Test.Tasty
import Test.Tasty.Golden
import Data.Vector.Generic
import Bio.Util.Bed
import Bio.Util.Overlap

main ∷ IO ()
main = defaultMain tests

tests :: TestTree
tests = testGroup "Tests" [unitTests]

unitTests ∷ TestTree
unitTests = testGroup "Unit tests"
  [ goldenVsString "Nonoverlapping features" "tests/test.txt" test1
  ]

parseTags ∷ [B.ByteString] → [Int]
parseTags xs = fmap (f . B.words) xs
    where
        f (_:s:e:_:_:strand:_) | strand == "+" = readInt s
                               | strand == "-" = readInt e 
                               | otherwise = readInt s
        f _ = -1

parseBed ∷ [B.ByteString] → [(Int, Int)]
parseBed xs = fmap (f . B.words) xs
    where
        f (_:s:e:_) = (readInt s, readInt e)
        f _ = (0,0)

test1 ∷ IO B.ByteString
test1 = do
    let f1 = "tests/f1.bed"
        f2 = "tests/f2.bed"
    feats ← B.readFile f1
    tags ← B.readFile f2
    let c = overlapFragment (parseBed $ B.lines feats) (parseBed $ B.lines tags)
    return $ B.unlines $ fmap (B.pack.show) $ toList c
