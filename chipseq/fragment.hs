{-# LANGUAGE OverloadedStrings, UnicodeSyntax, BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}

import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as V
import qualified Data.HashSet as S
import Data.Maybe
import Data.Bits
import Data.List
import System.Environment
import Graphics.Rendering.HPlot
import Data.Default
import Control.Concatenative
import Control.Parallel.Strategies

apprxCorr ∷ S.HashSet Int → S.HashSet Int → Int → Double
apprxCorr forwd rev d = fromIntegral $ S.foldl' f 0 rev
    where
        f ∷ Int → Int → Int
        f !acc r | any (`S.member` forwd) [r-d-4..r-d+4] = acc + 1
                 | otherwise = acc

type Interval = (Int, Int)

size ∷ Interval → Int
size (a,b) = b - a + 1

compare_ ∷ Int → Interval → Ordering
compare_ t (i, j) | t < i = LT
             | t > j = GT
             | otherwise = EQ

-- | return the overlaps of two lists of sorted non-overlapping intervals
overlap ∷ [Interval] → [Interval] → [Interval]
overlap ((h1, e1) : xs) ((h2, e2) : ys)
    -- case 1: h1 ------ e1
    --             h2 ------- e2
    | h2 > h1 && e1 >= h2 && e2 >= e1 = (h2, e1) : overlap xs ((e1,e2):ys)
    -- case 2:     -------
    --         -------
    | h1 > h2 && e2 >= h1 && e1 >= e2 = (h1, e2) : overlap ((e2,e1):xs) ys
    -- case 3:  -----
    --                ------
    | e1 < h2 = overlap xs ((h2,e2):ys)
    -- case 4:        ------
    --         ------
    | e2 < h1 = overlap ((h1,e1):xs) ys
    -- case 5:  ---------
    --            -----
    | h2 >= h1 && e2 <= e1 = (h2, e2) : overlap ((e2, e1):xs) ys
    -- case 6:    ----
    --          --------
    | h1 >= h2 && e1 <= e2 = (h1, e1) : overlap xs ((e1, e2):ys)
    | otherwise = error "overlap error"
overlap _ _ = []

overlap' ∷ [Interval] → [Interval] → [Interval]
overlap' a b = go a b []
    where
        go ((h1, e1) : xs) ((h2, e2) : ys) r
            | h2 > h1 && e1 >= h2 && e2 >= e1 = go xs ((e1,e2):ys) ((h2,e1):r)
            | h1 > h2 && e2 >= h1 && e1 >= e2 = go ((e2,e1):xs) ys ((h1,e2):r)
            | e1 < h2 = go xs ((h2,e2):ys) r
            | e2 < h1 = go ((h1,e1):xs) ys r
            | h2 >= h1 && e2 <= e1 = go ((e2, e1):xs) ys ((h2,e2):r)
            | h1 >= h2 && e1 <= e2 = go xs ((e1, e2):ys) ((h1,e1):r)
        go _ _ r = r

binarySearchWithBound ∷ G.Vector v Interval ⇒ (Int → Interval → Ordering) → Int → v Interval → Int → Int → Bool
binarySearchWithBound cmp t v i j
    | i > j = False
    | otherwise = let c = (i + j) `shiftR` 1
                  in case t `cmp` (v `G.unsafeIndex` c) of
        LT → binarySearchWithBound cmp t v i (c-1)
        GT → binarySearchWithBound cmp t v (c+1) j
        _ → True

binarySearch ∷ G.Vector v Interval ⇒ Int → v Interval → Bool
binarySearch t v = binarySearchWithBound compare_ t v 0 $ G.length v


toInt ∷ B.ByteString → Int
toInt = fst . fromJust . B.readInt

readBed ∷ [B.ByteString] → ([Int], [Int])
readBed = foldl' go ([],[])
    where
        go (!forwd, !rev) x = let (_:start:end:_:_:strand:_) = B.words x
                              in case strand of
                                  "+" → (toInt start : forwd, rev)
                                  "-" → (forwd, toInt end : rev)
                                  _ → error "read BED fail"


readInterval ∷ [B.ByteString] → [Interval]
readInterval = map f 
    where
        f x = let (_:start:end:_:_) = B.words x
              in (toInt start, toInt end)

countReads ∷ (Num a, G.Vector v Interval) ⇒ S.HashSet Int → v Interval → a
countReads rs m = S.foldl' f 0 rs
    where
        f acc i | binarySearch i m = acc + 1
                | otherwise = acc 

count_ ∷ (Num a, G.Vector v Interval) ⇒ S.HashSet Int → S.HashSet Int → v Interval → Int → a
count_ forwd rev m d = S.foldl' f 0 rev
    where
        f acc i = let i' = i - d
                  in if any (`S.member` forwd) [i'-4..i'+4] && binarySearch i' m
                        then acc + 1
                        else acc
    

corrAtD ∷ S.HashSet Int → S.HashSet Int → [Interval] → Int → Int → Double
corrAtD forwd rev m readlength d = (count / n - μ1 * μ2) / sqrt (v1 * v2)
    where
        dm = V.fromList $ overlap m $ map (both (\i → i + d - readlength + 1)) m
        n = fromIntegral $ V.foldl' (+) 0 $ V.map size dm
        μ1 = countReads forwd dm / n
        μ2 = countReads rev dm / n
        v1 = μ1 * (1 - μ1)
        v2 = μ2 * (1 - μ2)
        count = count_ forwd rev dm d

corrOverRange ∷ S.HashSet Int → S.HashSet Int → [Interval] → Int → [Int] → [Double]
corrOverRange forwd rev m readlength = parMap rseq (corrAtD forwd rev m readlength)

main ∷ IO ()
main = do
    [f1, f2] ← getArgs
    tags ← B.readFile f1
    maps ← B.readFile f2
    let (forwd, rev) = readBed $ B.lines tags
        intervals = readInterval $ tail $ B.lines maps
        forwd' = S.fromList forwd
        rev' = S.fromList rev
        readlength = 36
        x = [0,2..400]
        y = corrOverRange forwd' rev' intervals readlength x
        y' = parMap rseq (apprxCorr forwd' rev') x
--    _ ← plot' def [line def (Just $ map fromIntegral x, y)] "2.png"
    _ ← plot' def [line def (Just $ map fromIntegral x, y')] "1.png"
    return ()
