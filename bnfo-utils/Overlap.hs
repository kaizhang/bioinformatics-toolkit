{-# LANGUAGE OverloadedStrings, UnicodeSyntax, TemplateHaskell, BangPatterns #-}
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.IntervalMap.Strict as IM
import Data.List
import Data.Maybe
import System.Environment

countOverlap ∷ Ord a ⇒ 
               [(a, a)] -- ^ Ascending order list
             → [(a, a)]
             → [Int]
countOverlap xs ts = IM.elems $ foldl' search imap $ map toInterval ts
    where
        search m t = let intervals = fst $ unzip $ m `IM.intersecting` t
                     in foldl' (\ m0 x → IM.adjust (+1) x m0) m intervals
        imap = IM.fromAscList $ zip (map toInterval xs) [0,0..]
        toInterval (l, u) = IM.ClosedInterval l u

readInt ∷ B.ByteString → Int
readInt = fst.fromJust.B.readInt

readBed ∷ B.ByteString → [(Int, Int)]
readBed b = map (f.B.words) $ B.lines b
    where
        f (_:s:e:_) = (readInt s, readInt e)

main = do
    [f1, f2] ← getArgs
    feats ← B.readFile f1
    tags ← B.readFile f2
    mapM_ print $ countOverlap (readBed feats) $ readBed tags
