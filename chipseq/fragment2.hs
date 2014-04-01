{-# LANGUAGE OverloadedStrings, UnicodeSyntax, BangPatterns #-}

import qualified Data.ByteString.Lazy.Char8 as B
import Data.Maybe
import qualified Data.HashSet as S
import Data.List
import Data.Default
import Graphics.Rendering.HPlot

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

apprxCorr ∷ (S.HashSet Int, [Int]) → Int → Int
apprxCorr (forwd, rev) d = foldl' f 0 rev
    where
        f !acc r | any (`S.member` forwd) [r-d-10..r-d+10] = acc + 1
                 | otherwise = acc
    

main ∷ IO ()
main = do
    f ← B.getContents
    let (forwd, rev) = readBed $ B.lines f
        x = [0, 5..500]
        rev' = S.toList $ S.fromList rev
        y = map (fromIntegral.apprxCorr (S.fromList forwd, rev)) x

    plot' def [line def (Just $ map fromIntegral x, y)] "2.png"
    return ()

