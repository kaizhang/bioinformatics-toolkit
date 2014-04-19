{-# LANGUAGE OverloadedStrings, UnicodeSyntax, BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}

module ChIP (
      naiveCC
    , naiveCC'
    , corrOverRange
    , slidingAverage
    , readInterval
    , binarySearch
) where

import qualified Data.ByteString.Char8 as B
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as V
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import Data.Maybe
import Data.Bits
import Control.Concatenative
import Control.Parallel.Strategies
import Statistics.Sample (mean)

-- | fast relative cross-correlation with smoothing
apprxCorr ∷ S.HashSet Int → S.HashSet Int → Int → Int → Double
{-# INLINE apprxCorr #-}
apprxCorr forwd rev smooth d = fromIntegral $ S.foldl' f 0 rev
    where
        f ∷ Int → Int → Int
        f !acc r | any (`S.member` forwd) [r-d-smooth..r-d+smooth] = acc + 1
                 | otherwise = acc

naiveCC ∷ B.ByteString → [Int] → Double
naiveCC = naiveCCWithSmooth 0

naiveCC' ∷ B.ByteString → [Int] → Double
naiveCC' = naiveCCWithSmooth 4

naiveCCWithSmooth ∷ Int → B.ByteString → [Int] → Double
{-# INLINE naiveCCWithSmooth #-}
naiveCCWithSmooth smooth input range = mean.V.fromList.map (maxCC.snd) $ chrs
    where
        chrs = readBed input
        maxCC (forwd, rev) = fromIntegral.snd.maximum $ zip (parMap rpar (apprxCorr forwd rev smooth) range) range

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
{-# INLINE toInt #-}
toInt = fst . fromJust . B.readInt

readBed ∷ B.ByteString → [(B.ByteString, (S.HashSet Int, S.HashSet Int))]
{-# INLINE readBed #-}
readBed = map toSet.M.toList.M.fromListWith f.map parseLine.B.lines
    where
        parseLine ∷ B.ByteString → (B.ByteString, ([Int], [Int]))
        parseLine x = let (chr:start:end:_:_:strand:_) = B.words x
                      in case strand of
                          "+" → (chr, ([toInt start], []))
                          "-" → (chr, ([], [toInt end]))
                          _ → error "Unknown Strand!"
        f (a,b) (a',b') = (a++a', b++b')
        toSet (chr, (forwd, rev)) = (chr, (S.fromList forwd, S.fromList rev))


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
                  in if i' `S.member` forwd && binarySearch i' m
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

slidingAverage ∷ [Double] → Int → [Double]
slidingAverage xs windowSize = map f [0..n-1] 
    where
        xs' = V.fromList xs 
        n = V.length xs'
        f i = mean $ V.slice i' l xs'
            where
                i' = if i - windowSize >= 0 then i - windowSize else 0
                l = if i + windowSize <= n then windowSize else n - i

