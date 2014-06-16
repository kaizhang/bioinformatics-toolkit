module Bio.Seq 
    ( getSeqs
    , getSeq
    , getIndex
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict as M
import System.IO
import Data.List.Split
import Bio.Utils.Bed (readInt)

type IndexTable = M.HashMap B.ByteString Int
type OffSet = Int
type Query = (B.ByteString, Int, Int) -- zero based

getSeqs :: FilePath -> [Query] -> IO [B.ByteString]
getSeqs fl querys = do h <- openFile fl ReadMode
                       (index, offset) <- getIndex h
                       r <- mapM (getSeq h index offset) querys
                       hClose h
                       return r

getSeq :: Handle -> IndexTable -> OffSet -> Query -> IO B.ByteString
getSeq h index offset (chr, start, end) = do 
    hSeek h AbsoluteSeek (fromIntegral pos)
    B.hGet h (end - start + 1)
  where
    pos = offset + chrStart + start
    chrStart = M.lookupDefault (error "error") chr index

getIndex :: Handle -> IO (IndexTable, OffSet)
getIndex h = do header <- B.hGetLine h
                return (M.fromList.map f . chunksOf 2 . B.words $ header, B.length header + 1)
  where
    f [k, v] = (k, readInt v)
    f _ = error "error"
