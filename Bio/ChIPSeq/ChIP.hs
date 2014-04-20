{-# LANGUAGE OverloadedStrings, UnicodeSyntax, BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}

module Bio.ChIPSeq.ChIP (
      naiveCC
    , naiveCC'
) where

import qualified Data.ByteString.Char8 as B
import qualified Data.HashSet as S
import qualified Data.HashMap.Strict as M
import Data.List
import Bio.Util.Bed
import Control.Parallel.Strategies

readBed ∷ B.ByteString → [(B.ByteString, (S.HashSet Int, S.HashSet Int))]
{-# INLINE readBed #-}
readBed = map toSet.M.toList.M.fromListWith f.map parseLine.B.lines
    where
        parseLine ∷ B.ByteString → (B.ByteString, ([Int], [Int]))
        parseLine x = let (chr:start:end:_:_:strand:_) = B.words x
                      in case strand of
                          "+" → (chr, ([readInt start], []))
                          "-" → (chr, ([], [readInt end]))
                          _ → error "Unknown Strand!"
        f (a,b) (a',b') = (a++a', b++b')
        toSet (chr, (forwd, rev)) = (chr, (S.fromList forwd, S.fromList rev))

-- | fast relative cross-correlation with smoothing
apprxCorr ∷ S.HashSet Int → S.HashSet Int → Int → Int → Int
{-# INLINE apprxCorr #-}
apprxCorr forwd rev smooth d = S.foldl' f 0 rev
    where
        f ∷ Int → Int → Int
        f !acc r | any (`S.member` forwd) [r-d-smooth..r-d+smooth] = acc + 1
                 | otherwise = acc

naiveCC ∷ B.ByteString → [Int] → Int
naiveCC = naiveCCWithSmooth 0

naiveCC' ∷ B.ByteString → [Int] → Int
naiveCC' = naiveCCWithSmooth 4

naiveCCWithSmooth ∷ Int → B.ByteString → [Int] → Int
{-# INLINE naiveCCWithSmooth #-}
naiveCCWithSmooth smooth input range = maxCC.readBed $ input
    where
        cc ∷ [(B.ByteString, (S.HashSet Int, S.HashSet Int))] → Int → Int
        cc xs d = foldl' (+) 0 . map ((\(forwd, rev) → apprxCorr forwd rev smooth d).snd) $ xs
        maxCC xs = snd.maximum $ zip (parMap rpar (cc xs) range) range

