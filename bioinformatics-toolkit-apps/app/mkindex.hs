module Main where

import Bio.Seq.IO (mkIndex)
import System.Environment

main :: IO ()
main = do
    (outF: inputs)  <- getArgs
    mkIndex inputs outF
