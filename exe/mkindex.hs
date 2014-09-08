module Main where

import Bio.Seq.Query (mkIndex)
import System.Environment

main :: IO ()
main = do [inF, outF] <- getArgs
          mkIndex inF outF

