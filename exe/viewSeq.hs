module Main where

import qualified Data.ByteString.Char8 as B
import Bio.Seq
import Bio.Seq.IO
import System.Environment

main :: IO ()
main = do
    [fl, chr, start, end] <- getArgs
    g <- pack fl
    s <- getSeqs g [(B.pack chr, read start, read end)] :: IO [Either String (DNA IUPAC)]
    case head s of
        Left err -> error err
        Right x -> B.putStrLn . toBS $ x
