{-# LANGUAGE OverloadedStrings #-}

import Criterion.Main
import System.Random
import qualified Data.ByteString.Char8 as B
import Bio.Motif
import Bio.Seq
import Data.Default.Generics

dna :: DNA Basic
dna = fromBS $ B.pack $ map f $ take 5000 $ randomRs (0, 3) (mkStdGen 2)
  where
    f :: Int -> Char
    f x = case x of
        0 -> 'A'
        1 -> 'C'
        2 -> 'G'
        3 -> 'T'
        _ -> undefined
       

pwm :: PWM
pwm = readPWM "0.3 0.3 0.3 0.1\n0 0.5 0 0.5\n0.1 0.2 0.5 0.3\n0.1 0.1 0.1 0.7\n0 0 0 1\n0.25 0.25 0.25 0.25"

main :: IO ()
main = defaultMain [
      bench "motif score" $ nf (scores def pwm) dna 
    ]
