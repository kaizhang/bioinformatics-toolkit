{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

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

import Control.Exception (bracket)
import Conduit
import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Lazy as M
import Data.List.Split
import System.IO

import Bio.Seq
import Bio.Utils.Misc (readInt)
import Bio.Data.Fasta (fastaReader)

-- | The first 2048 bytes are header. Header consists of a magic string,
-- followed by chromosome information. Example:
-- <HASKELLBIOINFORMATICS>\nCHR1 START SIZE
data Genome = Genome !Handle !IndexTable !Int

type IndexTable = M.HashMap B.ByteString (Int, Int)

magic :: B.ByteString
magic = "<HASKELLBIOINFORMATICS_7d2c5gxhg934>"

openGenome :: FilePath -> IO Genome
openGenome fl = do
    h <- openFile fl ReadMode
    sig <- B.hGetLine h
    if sig == magic
        then do
            header <- B.hGetLine h
            return $ Genome h (getIndex header) (B.length sig + B.length header + 2)
        else error "Bio.Seq.Query.openGenome: Incorrect format"
  where
    getIndex = M.fromList . map f . chunksOf 3 . B.words
      where
        f [k, v, l] = (k, (readInt v, readInt l))
        f _ = error "error"
{-# INLINE openGenome #-}

closeGenome :: Genome -> IO ()
closeGenome (Genome h _ _) = hClose h
{-# INLINE closeGenome #-}

withGenome :: FilePath -> (Genome -> IO a) -> IO a
withGenome fl fn = bracket (openGenome fl) closeGenome fn
{-# INLINE withGenome #-}

-- | A query is represented by a tuple: (chr, start, end) and is
-- zero-based index, half-close-half-open
type Query = (B.ByteString, Int, Int)

-- | Retrieve sequence.
getSeq :: BioSeq s a => Genome -> Query -> IO (Either String (s a))
getSeq (Genome h index headerSize) (chr, start, end) = case M.lookup chr index of
    Just (chrStart, chrSize) ->
        if end > chrSize
            then return $ Left $ "Bio.Seq.getSeq: out of index: " ++ show end ++
                ">" ++ show chrSize
            else do
                hSeek h AbsoluteSeek $ fromIntegral $ headerSize + chrStart + start
                fromBS <$> B.hGet h (end - start)
    _ -> return $ Left $ "Bio.Seq.getSeq: Cannot find " ++ show chr
{-# INLINE getSeq #-}

-- | Retrieve whole chromosome.
getChrom :: Genome -> B.ByteString -> IO (Either String (DNA IUPAC))
getChrom g chr = case lookup chr chrSize of
    Just s -> getSeq g (chr, 0, s)
    _ -> return $ Left "Unknown chromosome"
  where
    chrSize = getChrSizes g
{-# INLINE getChrom #-}

-- | Retrieve chromosome size information.
getChrSizes :: Genome -> [(B.ByteString, Int)]
getChrSizes (Genome _ table _) = map (\(k, (_, l)) -> (k, l)) . M.toList $ table
{-# INLINE getChrSizes #-}

-- | Indexing a genome.
mkIndex :: [FilePath]    -- ^ A list of fasta files. Each file can have multiple
                         -- chromosomes.
        -> FilePath      -- ^ output file
        -> IO ()
mkIndex fls outFl = withFile outFl WriteMode $ \outH -> do
    (chrs, dnas) <- (unzip . concat) <$> mapM readSeq fls
    B.hPutStr outH $ B.unlines [magic, mkHeader chrs]
    mapM_ (B.hPutStr outH) dnas
  where
    readSeq fl = runResourceT $ fastaReader fl =$= conduit $$ sinkList
      where
        conduit = awaitForever $ \(chrName, seqs) -> do
            let dna = B.concat seqs
            yield ((head $ B.words chrName, B.length dna), dna)
{-# INLINE mkIndex #-}

mkHeader :: [(B.ByteString, Int)] -> B.ByteString
mkHeader xs = B.unwords.fst $ foldl f ([], 0) xs
  where
    f (s, i) (s', i') = (s ++ [s', B.pack $ show i, B.pack $ show i'], i + i')
{-# INLINE mkHeader #-}
