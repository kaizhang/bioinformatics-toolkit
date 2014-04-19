{-# LANGUAGE OverloadedStrings, UnicodeSyntax, TemplateHaskell, BangPatterns #-}

module Bio.Util.Overlap (
      overlapFragment
    , overlapNucl
) where

import qualified Data.ByteString.Char8 as B
import qualified Data.IntervalMap.Strict as IM
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as VM
import Bio.Util.Bed
import Control.Monad

overlapFragment, overlapNucl ∷ 
                  [(Int, Int)] -- ^ Ascending order list
                → [(Int, Int)] -- ^ tags in any order
                → V.Vector Int
overlapFragment xs ts = V.create (VM.replicate n 0 >>= go ts)
    where
        n = length xs
        iMap = IM.fromAscList $ zip (map toInterval xs) [0..]
        go ts' v = do
            forM_ ts' (\x → do
                let indices = snd.unzip.IM.intersecting iMap.toInterval $ x
                forM_ indices (\i → VM.write v i . (+1) =<< VM.read v i)
                )
            return v

overlapNucl xs ts = V.create (VM.replicate n 0 >>= go ts)
    where
        n = length xs
        iMap = IM.fromAscList $ zip (map toInterval xs) [0..]
        go ts' v = do
            forM_ ts' (\x → do
                let intervals = IM.intersecting iMap.toInterval $ x
                forM_ intervals (\interval → do 
                    let i = snd interval
                        nucl = overlap x . fst $ interval
                    VM.write v i . (+nucl) =<< VM.read v i
                    )
                )
            return v
        overlap (l, u) (IM.ClosedInterval l' u') 
            | l' >= l = if u' <= u then u'-l'+1 else u-l'+1
            | otherwise = if u' <= u then u'-l+1 else u-l+1
        overlap _ _ = 0

toInterval ∷ (a, a) → IM.Interval a
{-# INLINE toInterval #-}
toInterval (l, u) = IM.ClosedInterval l u

readBed ∷ B.ByteString → [(Int, Int)]
readBed b = map (f.B.words) $ B.lines b
    where
        f (_:s:e:_) = (readInt s, readInt e)
