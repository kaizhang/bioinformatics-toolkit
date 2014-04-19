{-# LANGUAGE OverloadedStrings, UnicodeSyntax #-}

import qualified Data.ByteString.Lazy.Char8 as B
import System.Environment
import Graphics.Rendering.HPlot
import Data.Default
import Control.Parallel.Strategies
import qualified Data.HashSet as S
import ChIP
import Control.Lens hiding (both)
import Control.Concatenative
import Spline

run ∷ [FilePath] → IO ()
run [f1,f2,f3] = do
    tags ← B.readFile f1
    maps ← B.readFile f2
    let (forwd, rev) = both S.fromList $ readBed $ B.lines tags
        m = readInterval $ tail $ B.lines maps
        x = [0,1..400]
        y = corrOverRange forwd rev m (read f3) x
        y' = slidingAverage y 15
    _ ← plot' def [
        line def (Just $ map fromIntegral x, y)
--        line (opacity .~ 0.7 $ col .~ "green" $ def) (Just $ map fromIntegral x, y')
        ] "1.png"
    return ()

run [f1,f2] = do
    tags ← B.readFile f1
    let (forwd, rev) = both S.fromList $ readBed $ B.lines tags
        x = [0,1..400]
        y = parMap rseq (apprxCorr forwd rev (read f2)) x
    (x',y') ← smoothSpline (zip (map fromIntegral x) y) 0.6
    _ ← plot' def [
        line def (Just $ map fromIntegral x, y)
--        line (col .~ "black" $ def) (Just x', y')
        ] "1.png"
    return ()

run _ = error "incorrect arguments"

main ∷ IO ()
main = do
    args ← getArgs
    run args
