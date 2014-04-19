{-# LANGUAGE OverloadedStrings, UnicodeSyntax, BangPatterns #-}

module Bio.Util.Bed (
      readInt
    , readDouble
    , bins
    , binBySize
) where

import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.ByteString.Lex.Double as L
import qualified Data.ByteString.Lex.Lazy.Double as LL
import Data.Maybe

class ToNum a where
    readInt ∷ a → Int
    readDouble ∷ a → Double

instance ToNum B.ByteString where
    readInt = fst.fromJust.B.readInt
    readDouble = fst . fromJust . L.readDouble

instance ToNum BL.ByteString where
    readInt = fst.fromJust.BL.readInt
    readDouble = fst . fromJust . LL.readDouble

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
