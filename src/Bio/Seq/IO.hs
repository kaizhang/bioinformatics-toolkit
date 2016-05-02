{-# LANGUAGE BangPatterns #-}

module Bio.Seq.IO
    ( Genome
    , openGenome
    , closeGenome
    , withGenome
    , getSeq
    , getChrom
    , getChrSizes
    , mkIndex
    ) where

import Bio.Seq
import Bio.Utils.Misc (readInt)
import Control.Exception (bracket)
import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Lazy as M
import Data.List.Split
import System.IO

-- | The first 2048 bytes are header. Header consists of a magic string,
-- followed by chromosome information. Example:
-- <HASKELLBIOINFORMATICS>\nCHR1 START SIZE
data Genome = Genome !Handle !IndexTable

type IndexTable = M.HashMap B.ByteString (Int, Int)

headerSize :: Int
headerSize = 2048

magic :: String
magic = "<HASKELLBIOINFORMATICS>"

openGenome :: FilePath -> IO Genome
openGenome fl = do
    h <- openFile fl ReadMode
    sig <- hGetLine h
    if sig == magic
        then Genome h <$> readIndex h
        else error "Bio.Seq.Query.openGenome: Incorrect format"
  where
    readIndex h = do
        header <- B.hGetLine h
        return $ M.fromList . map f . chunksOf 3 . B.words $ header
      where
        f [k, v, l] = (k, (readInt v, readInt l))
        f _ = error "error"
{-# INLINE openGenome #-}

closeGenome :: Genome -> IO ()
closeGenome (Genome h _) = hClose h
{-# INLINE closeGenome #-}

withGenome :: FilePath -> (Genome -> IO a) -> IO a
withGenome fl fn = bracket (openGenome fl) closeGenome fn
{-# INLINE withGenome #-}

type Query = (B.ByteString, Int, Int) -- (chr, start, end), zero-based index, half-close-half-open

getSeq :: BioSeq s a => Genome -> Query -> IO (Either String (s a))
getSeq (Genome h index) (chr, start, end) = case M.lookup chr index of
    Just (chrStart, chrSize) ->
        if end > chrSize
            then return $ Left $ "Bio.Seq.getSeq: out of index: " ++ show end ++
                ">" ++ show chrSize
            else do
                hSeek h AbsoluteSeek $ fromIntegral $ headerSize + chrStart + start
                (Right . fromBS) <$> B.hGet h (end - start)
    _ -> return $ Left $ "Bio.Seq.getSeq: Cannot find " ++ show chr
{-# INLINE getSeq #-}

getChrom :: Genome -> B.ByteString -> IO (Either String (DNA IUPAC))
getChrom g chr = case lookup chr chrSize of
    Just s -> getSeq g (chr, 0, s)
    _ -> return $ Left "Unknown chromosome"
  where
    chrSize = getChrSizes g
{-# INLINE getChrom #-}

getChrSizes :: Genome -> [(B.ByteString, Int)]
getChrSizes (Genome _ table) = map (\(k, (_, l)) -> (k, l)) . M.toList $ table
{-# INLINE getChrSizes #-}

-- | indexing a genome.
mkIndex :: [FilePath]    -- ^ fasta files representing individual chromosomes
        -> FilePath      -- ^ output file
        -> IO ()
mkIndex fls outFl = withFile outFl WriteMode $ \outH -> do
    hPutStr outH $ magic ++ "\n" ++ replicate 2024 '#'  -- header size: 1024
    chrs <- mapM (write outH) fls
    hSeek outH AbsoluteSeek 24
    B.hPutStrLn outH $ mkHeader chrs
  where
    write handle fl = withFile fl ReadMode $ \h -> do
        fastaHeader <- B.hGetLine h
        n <- loop 0 h
        return (B.tail fastaHeader, n)
      where
        loop !n h' = do eof <- hIsEOF h'
                        if eof then return n
                               else do l <- B.hGetLine h'
                                       B.hPutStr handle l
                                       loop (n + B.length l) h'
{-# INLINE mkIndex #-}

mkHeader :: [(B.ByteString, Int)] -> B.ByteString
mkHeader xs = B.unwords.fst $ foldl f ([], 0) xs
  where
    f (s, i) (s', i') = (s ++ [s', B.pack $ show i, B.pack $ show i'], i + i')
{-# INLINE mkHeader #-}
