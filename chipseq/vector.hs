{-# LANGUAGE OverloadedStrings, UnicodeSyntax, TemplateHaskell, BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}

import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.HashSet as S
import Control.Monad
import Data.Maybe
import Data.List
import System.Environment
import Graphics.Rendering.HPlot
import Data.Default

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


readInterval ∷ [B.ByteString] → [(Int, Int)]
readInterval = map f 
    where
        f x = let (chr:start:end:m:_) = B.words x
              in (toInt start, toInt end)

mkMap ∷ G.Vector v Bool ⇒ [(Int, Int)] → Int → v Bool
mkMap intervals chrSize = G.create (GM.replicate (chrSize+1) False >>= go intervals)
    where
        go (x:xs) v = do
            forM_ [fst x .. snd x] (\i → 
                GM.write v i True)
            go xs v
        go _ v = return v

-- | Is given position mappable at reverse strand shifted by d
mappableR ∷ G.Vector v Bool ⇒ v Bool → Int → Int → Int → Bool
mappableR m readlength i d = let shiftedPos = i + d - readlength + 1
                                 l = G.length m
                               in if shiftedPos >= 0 && shiftedPos < l
                                     then G.unsafeIndex m shiftedPos
                                     else False

doubleMappable ∷ G.Vector v Bool ⇒ v Bool → Int → Int → Int → Bool
doubleMappable m readlength i d | i < G.length m = G.unsafeIndex m i && mappableR m readlength i d
                                | otherwise = False

-- | number of double mappable regions given d
numOfDM ∷ (Num a, G.Vector v Bool) ⇒ v Bool → Int → Int → a
numOfDM m readlength d = foldl' f 0 [0..l]
    where
        l = G.length m
        f acc i | doubleMappable m readlength i d = acc + 1
                | otherwise = acc

countReadsDM ∷ (Num a, G.Vector v Bool) ⇒ S.HashSet Int → v Bool → Int → Int → a
countReadsDM rs m readlength d = S.foldl' f 0 rs
    where
        f acc i | doubleMappable m readlength i d = acc + 1
                | otherwise = acc 

count_ ∷ (Num a, G.Vector v Bool) ⇒ S.HashSet Int → S.HashSet Int → v Bool → Int → Int → a
count_ forwd rev m readlength d = S.foldl' f 0 rev
    where
        f acc i = let i' = i - d
                  in if (i' `S.member` forwd && doubleMappable m readlength i' d)
                        then acc + 1
                        else acc
    

corrAtD ∷ G.Vector v Bool ⇒ S.HashSet Int → S.HashSet Int → v Bool → Int → Int → Double
corrAtD forwd rev m readlength d = (count / n - μ1 * μ2) / sqrt (v1 * v2)
    where
        n = numOfDM m readlength d
        μ1 = countReadsDM forwd m readlength d / n
        μ2 = countReadsDM rev m readlength d / n
        v1 = μ1 * (1 - μ1)
        v2 = μ2 * (1 - μ2)
        count = count_ forwd rev m readlength d

corrOverRange ∷ G.Vector v Bool ⇒ S.HashSet Int → S.HashSet Int → v Bool → Int → [Int] → [Double]
corrOverRange forwd rev m readlength = map (corrAtD forwd rev m readlength)

main ∷ IO ()
main = do
    [f1, f2] ← getArgs
    reads ← B.readFile f1
    maps ← B.readFile f2
    let (forwd, rev) = readBed $ B.lines reads
        intervals = readInterval $ tail $ B.lines maps
        chrSize = 166650296
        forwd' = S.fromList forwd
        rev' = S.fromList rev
        m = mkMap intervals chrSize ∷ V.Vector Bool
        readlength = 36
        x = [1,3..200]
        y = corrOverRange forwd' rev' m readlength x
    plot' def [line def (Just $ map fromIntegral x, y)] "2.png"
    return ()
