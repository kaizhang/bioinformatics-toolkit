module Main where

import Bio.Seq.IO (mkIndex)
import Shelly
import qualified Data.Text as T
import System.Environment

main :: IO ()
main = do [inDir, outF] <- getArgs
          fls <- shelly . lsT . fromText . T.pack $ inDir
          mkIndex (map T.unpack fls) outF
