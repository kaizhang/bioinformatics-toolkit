{-# LANGUAGE OverloadedStrings, UnicodeSyntax #-}

module Spline (
    smoothSpline
) where

import Shelly
import qualified Data.Text as T
import Data.List

smoothSpline ∷ [(Double, Double)] → Double → IO ([Double], [Double])
smoothSpline datas i = procResult =<< shelly (silently $ run "Rscript" ["-e", script])
    where
        script = T.intercalate "" ["a <- smooth.spline(", x, ", ", y, ", spar=", λ, ");cat(a$x);cat('####');cat(a$y)"]
        (x', y') = unzip datas
        x = T.pack $ "c(" ++ intercalate ", " (map show x') ++ ")"
        y = T.pack $ "c(" ++ intercalate ", " (map show y') ++ ")"
        λ = T.pack $ show i
        procResult r = (return . f . map (map (read . T.unpack) . tail . T.words) . T.splitOn "####") r
            where
                f [a, b] = (a,b)
                f _ = error $ T.unpack r
