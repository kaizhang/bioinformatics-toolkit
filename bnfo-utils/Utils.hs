{-# LANGUAGE OverloadedStrings, UnicodeSyntax, TemplateHaskell, BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Rank2Types #-}

module BNFO.Utils (
      (+>)
    , readInt
    , readDouble
    , getFields
    , bins
    , binBySize
    , setField
    , toMap
    , binarySearchLE
    , binarySearchGE
    , countOverlap
    , countNonOverlap
    , countCNN
) where

import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lex.Double as L
import qualified Data.HashMap.Strict as M
import qualified Data.IntervalMap.Strict as IM
import Data.Maybe
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import Control.Monad
import Control.Concatenative
import Data.Bits


-- zero-index
getField ∷ B.ByteString → Int → B.ByteString
getField s n = (B.split '\t' s)!!n

(+>) = getField

getFields ∷ B.ByteString → [Int] → [B.ByteString]
getFields gff fs = fmap ((B.split '\t' gff)!!) fs

toMap ∷ B.ByteString → Int → [Int] → M.HashMap B.ByteString [B.ByteString]
toMap input k fields = M.fromList $ map f $ B.lines input
    where f line = let items = B.split '\t' line
                   in ((items!!(k-1)), map (\ x → items!!(x-1)) fields)

-- zero-index
setField ∷ Int → B.ByteString → B.ByteString → B.ByteString
setField n field str = let fields = B.split '\t' str
                       in B.intercalate "\t" $ take n fields ++ (field : drop (n+1) fields)

readInt ∷ B.ByteString → Int
readInt = fst.fromJust.B.readInt

readDouble ∷ B.ByteString → Double
readDouble = fst . fromJust . L.readDouble

-- | divide a given region into fixed size fragments
binBySize ∷ (Int, Int) → Int → [(Int, Int)]
binBySize (start, end) step =
    let binNum = ceiling $ fromIntegral (end - start + 1) / fromIntegral step
    in take binNum $ zip [start,start+step..] [start+step-1,start+2*step-1..]

-- | divide a given region into k equal size sub-regions
bins ∷ (Int, Int) → Int → [(Int, Int)]
bins (start, end) binNum = 
    let step = ceiling $ fromIntegral (end - start + 1) / fromIntegral binNum
    in take binNum $ zip [start,start+step..] [start+step-1,start+2*step-1..]

-- | read count with non-overlapping features
countNonOverlap ∷ (Ord a, G.Vector v Int) ⇒
                  [(a, a)]  -- ^ SORTED starting and ending locations of features
                → [a]         -- ^ starting locations of tags/reads
                → v Int
countNonOverlap feats tags = count_ feats tags nonOverlapHelp

-- | read count with overlapping features
countOverlap ∷ (Ord a, G.Vector v Int) ⇒
               [(a, a)]
             → [a]
             → v Int
countOverlap feats tags = count_ feats tags overlapHelp

-- | count reads on continuous non-overlapping features
countCNN ∷ (G.Vector v Int) ⇒
           Int  -- ^ Chromosome size
         → Int  -- ^ Bin size
         → [Int] -- ^ tags/reads
         → v Int
{-# INLINE countCNN #-}
countCNN chrSize binSize tags = G.create (GM.replicate binNum 0 >>= go tags)
    where
        binNum = ceiling $ fromIntegral chrSize / fromIntegral binSize
        go xs bin = do
            forM_ xs $ (\x → do
                let i = (x - 1) `div` binSize
                GM.unsafeWrite bin i . (+1) =<< GM.read bin i)
            return bin
        

-- | general read count function
count_ ∷ (Ord a, G.Vector v Int) ⇒
               [(a, a)]                  -- ^ SORTED starting and ending locations of features
             → [a]                       -- ^ starting locations of tags/reads
             → ([(a, a)] → a → [Int])    -- ^ core function
             → v Int
{-# INLINE count_ #-}
count_ feats tags locate = G.create (GM.replicate n 0 >>= go tags)
            where
                n = length tags
                go tags' xs = do
                    mapM_ count tags'
                    return xs
                        where
                            count = mapM_ 
                                (\i → GM.unsafeWrite xs i . (+1) =<< GM.unsafeRead xs i ).locate feats

nonOverlapHelp ∷ Ord a ⇒ [(a, a)] → a → [Int]
{-# INLINE nonOverlapHelp #-}
nonOverlapHelp feats t = fromMaybe [] $ do 
    i ← binarySearchLE starts t
    guard (G.unsafeIndex ends i >= t)
    return [i]
        where
            (starts, ends) = both V.fromList $ unzip feats

overlapHelp ∷ Ord a ⇒ [(a, a)] → a → [Int]
{-# INLINE overlapHelp #-}
overlapHelp feats t = map snd $ m `IM.containing` t
    where
        m = IM.fromAscList $ zip (map toInterval feats) [0..]
        toInterval (l, u) = IM.ClosedInterval l u

-- | return the index of the largest element that is less than or equal to given value
binarySearchLE ∷ (Ord a, G.Vector v a) ⇒ v a → a → Maybe Int
binarySearchLE xs b = binarySearchLByBound compare xs 0 (G.length xs - 1) b

-- | return the index of the largest element that is greater than or equal to given value
binarySearchGE ∷ (Ord a, G.Vector v a) ⇒ v a → a → Maybe Int
binarySearchGE xs b = binarySearchRByBound compare xs 0 (G.length xs - 1) b

binarySearchLByBound ∷ G.Vector v a ⇒ (a → a → Ordering) → v a → Int → Int → a → Maybe Int
binarySearchLByBound cmp xs l u b = toMaybe $ go l u
    where
        go i j | i >= j = case cmp (xs `G.unsafeIndex` j) b of
                                 GT → j - 1
                                 _ → j
               | otherwise = let c = (i + j) `shiftR` 1
                                 e = xs `G.unsafeIndex` c
                             in case cmp e b of
                                 LT → go (c+1) j
                                 GT → go i (c-1)
                                 _ → go i c
        toMaybe x | x < l = Nothing
                  | otherwise = Just x
            
binarySearchRByBound ∷ G.Vector v a ⇒ (a → a → Ordering) → v a → Int → Int → a → Maybe Int
binarySearchRByBound cmp xs l u b = toMaybe $ go l u
    where
        go i j | i >= j = case cmp (xs `G.unsafeIndex` j) b of
                                 LT → j + 1
                                 _ → j
               | otherwise = let c = (i + j + 1) `shiftR` 1
                                 e = xs `G.unsafeIndex` c
                             in case cmp e b of
                                 LT → go (c+1) j
                                 GT → go i (c-1)
                                 _ → go c j
        toMaybe x | x > u = Nothing
                  | otherwise = Just x
            
