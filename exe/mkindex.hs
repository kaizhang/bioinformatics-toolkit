module Main where

import Bio.Seq.IO (mkIndex)
import System.Environment

main :: IO ()
main = do [inF, outF] <- getArgs
          mkIndex inF outF

