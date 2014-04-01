{-# LANGUAGE OverloadedStrings, UnicodeSyntax, TemplateHaskell, BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}

import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Generic.Mutable as GM
import Control.Monad
import Data.Maybe
import Data.List
import System.Environment

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

storeTags ∷ G.Vector v Bool ⇒ [Int] → Int → v Bool
storeTags tags chrSize = G.create (GM.replicate chrSize False >>= go tags)
    where
        go (x:xs) v = do
            GM.write v x True
            go xs v
        go _ v = return v

mkMap ∷ G.Vector v Bool ⇒ [(Int, Int)] → Int → v Bool
mkMap intervals chrSize = G.create (GM.replicate chrSize False >>= go intervals)
    where
        go (x:xs) v = do
            forM_ [fst x .. snd x] (\i → 
                GM.write v i True)
            go xs v
        go _ v = return v

mkDMap ∷ G.Vector v Bool ⇒ v Bool → Int → Int → v Bool
mkDMap m d readlength = G.create (GM.replicate (G.length m) False >>= go)
    where
        go v = do
            forM_ [0..G.length m - 1] (\i → do
                let shift = i + d -readlength + 1
                when ((shift >= 0) && G.unsafeIndex m i && G.unsafeIndex m shift) $
                    GM.write v i True)
            return v

intersection ∷ G.Vector v Bool ⇒ v Bool → v Bool → v Bool
intersection a b = G.zipWith (&&) a b

numOfOnes ∷ (Num a, G.Vector v Bool) ⇒ v Bool → a
numOfOnes = fromIntegral . G.length . G.filter id

corr ∷ G.Vector v Bool ⇒ v Bool → v Bool → v Bool → Double
corr forwd rev m = (count / n - μ1 * μ2) / sqrt (v1 * v2)
    where
        n = numOfOnes m
        forwd' = intersection forwd m
        rev' = intersection rev m
        μ1 = numOfOnes forwd' / n
        μ2 = numOfOnes rev' / n
        v1 = μ1 * (1 - μ1)
        v2 = μ2 * (1 - μ2)
        count = numOfOnes $ intersection forwd' rev'

corrAtD ∷ G.Vector v Bool ⇒ v Bool → v Bool → v Bool → Int → Int → Double
corrAtD forwd rev m readlength d = result `seq` result
    where
        result = corr forwd rev dm
        dm = mkDMap m d readlength

corrOverRange ∷ G.Vector v Bool ⇒ v Bool → v Bool → v Bool → Int → [Int] → [Double]
corrOverRange forwd rev m readlength = map (corrAtD forwd rev m readlength)

main ∷ IO ()
main = do
    [f1, f2] ← getArgs
    reads ← B.readFile f1
    maps ← B.readFile f2
    let (forwd, rev) = readBed $ B.lines reads
        intervals = readInterval $ tail $ B.lines maps
        chrSize = 166650296
        forwd' = storeTags forwd chrSize ∷ V.Vector Bool
        rev' = storeTags rev chrSize ∷ V.Vector Bool
        m = mkMap intervals chrSize ∷ V.Vector Bool
        readlength = 36
        drange = [1..1]
    print $ corrOverRange forwd' rev' m readlength drange
