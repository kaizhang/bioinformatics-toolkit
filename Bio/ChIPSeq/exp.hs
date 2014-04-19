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
import Data.List

run ∷ [FilePath] → IO ()
run [f1] = do
    tags ← B.readFile f1
    let (forwd, rev) = both S.fromList $ readBed $ B.lines tags
        x = [0,1..400]
        ys = map (\i → parMap rseq (apprxCorr forwd rev i) x) [0,2..20]
--        y = map sum $ transpose $ zipWith (zipWith (/)) (tail ys) ys
--    (x',y') ← smoothSpline (zip (map fromIntegral x) y) 0.6
    _ ← plot' def (
        map (\y → line (opacity .~ 0.8 $ def) (Just $ map fromIntegral x, y)) ys
--        line (opacity .~ 0.8 $ def) (Just $ map fromIntegral x, y)
--        line (col .~ "black" $ def) (Just x', y')
        ) "2.png"
    return ()

run _ = error "incorrect arguments"

main ∷ IO ()
main = do
    args ← getArgs
    run args
