{-# LANGUAGE BangPatterns #-}

module Bio.Seq.Query
    ( Genome
    , GenomeH
    , gHOpen
    , gHClose
    , pack
    , getSeqs
    , getSeq
    , readIndex
    , mkIndex
    ) where

import Bio.Seq
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import Shelly hiding (FilePath)
import qualified Data.HashMap.Lazy as M
import System.IO
import Data.List.Split
import Bio.Utils.Misc (readInt)
import Control.Monad

-- | The first 2048 bytes are header. Header consists of a fingerprint string,
-- followed by chromosome information. Example:
-- <HASKELLBIOINFORMATICS>\nCHR1 START SIZE
newtype Genome = G FilePath

newtype GenomeH = GH Handle

fingerprint :: String
fingerprint = "<HASKELLBIOINFORMATICS>"

pack :: FilePath -> IO Genome
pack fl = withFile fl ReadMode f
  where f h = do l <- hGetLine h
                 if l == fingerprint
                    then return $ G fl
                    else error "Bio.Seq.Query.pack: Incorrect format"

gHOpen :: Genome -> IO GenomeH
gHOpen (G fl) = do h <- openFile fl ReadMode
                   return $ GH h

gHClose :: GenomeH -> IO ()
gHClose (GH h) = hClose h

type IndexTable = M.HashMap B.ByteString (Int, Int)

type Query = (B.ByteString, Int, Int) -- (chr, start, end), zero based

getSeqs :: BioSeq s a => Genome -> [Query] -> IO [s a]
getSeqs g querys = do gH <- gHOpen g
                      index <- readIndex gH
                      r <- mapM (getSeq gH index) querys
                      gHClose gH
                      return r

getSeq :: BioSeq s a => GenomeH -> IndexTable -> Query -> IO (s a)
getSeq (GH h) index (chr, start, end) = do 
    hSeek h AbsoluteSeek (fromIntegral pos)
    liftM fromBS $ B.hGet h (end - start + 1)
  where
    pos = headerSize + chrStart + start
    chrStart = fst $ M.lookupDefault (error $ "Bio.Seq.getSeq: Cannot find " ++ show chr) chr index
    headerSize = 2048

readIndex :: GenomeH -> IO IndexTable
readIndex (GH h) = do header <- B.hGetLine h >> B.hGetLine h
                      return $ M.fromList . map f . chunksOf 3 . B.words $ header
  where
    f [k, v, l] = (k, (readInt v, readInt l))
    f _ = error "error"

-- | indexing a genome.
mkIndex :: FilePath    -- ^ a directory containing fasta files for each chromosome
        -> FilePath    -- ^ output file
        -> IO ()
mkIndex dir outFl = do 
    outH <- openFile outFl WriteMode
    hPutStr outH $ fingerprint ++ "\n" ++ replicate 2024 '#'  -- header size: 1024
    fls <- (fmap.map) T.unpack $ shelly $ lsT (fromText $ T.pack dir)
    chrs <- mapM (write outH) fls
    hSeek outH AbsoluteSeek 24
    B.hPutStrLn outH $ mkHeader chrs
    hClose outH                      
  where
    write handle fl = do h <- openFile fl ReadMode
                         fastaHeader <- B.hGetLine h
                         n <- loop 0 h
                         hClose h
                         return (B.tail fastaHeader, n)
      where
        loop !n h' = do eof <- hIsEOF h'
                        if eof then return n
                               else do l <- B.hGetLine h'
                                       B.hPutStr handle l
                                       loop (n + B.length l) h'

mkHeader :: [(B.ByteString, Int)] -> B.ByteString
mkHeader xs = B.unwords.fst $ foldl f ([], 0) xs
  where
    f (s, i) (s', i') = (s ++ [s', B.pack $ show i, B.pack $ show i'], i + i')
{-# INLINE mkHeader #-}
