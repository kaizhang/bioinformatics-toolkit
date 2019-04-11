{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedStrings #-}
module Bio.Data.Fastq
    ( Fastq(..)
    , streamFastqGzip
    , streamFastq
    , sinkFastqGzip
    , sinkFastq
    , parseFastqC
    , fastqToByteString
    , qualitySummary
    , trimPolyA
    ) where

import           Conduit
import Data.Conduit.Zlib (ungzip, multiple, gzip)
import           Control.Monad         (when)
import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy as BL
import qualified Data.ByteString as BS
import           Data.Maybe            (isJust)
import qualified Data.Attoparsec.ByteString as A
import Data.Attoparsec.ByteString.Char8
import Data.Conduit.Attoparsec

-- | A FASTQ file normally uses four lines per sequence.
--
--     * Line 1 begins with a '@' character and is followed by a sequence
--       identifier and an optional description (like a FASTA title line).
--
--     * Line 2 is the raw sequence letters.
--
--     * Line 3 begins with a '+' character and is optionally followed by the
--       same sequence identifier (and any description) again.
--
--     * Line 4 encodes the quality values for the sequence in Line 2, and must
--       contain the same number of symbols as letters in the sequence.
data Fastq = Fastq
    { fastqSeqId   :: B.ByteString
    , fastqSeq     :: B.ByteString
    , fastqSeqQual :: B.ByteString
    } deriving (Show, Eq)

-- | Read gzipped fastq file.
streamFastqGzip :: (PrimMonad m, MonadThrow m, MonadResource m) 
                => FilePath -> ConduitT i Fastq m ()
streamFastqGzip fl = sourceFileBS fl .| multiple ungzip .| parseFastqC

streamFastq :: (MonadResource m, MonadThrow m)
            => FilePath -> ConduitT i Fastq m ()
streamFastq fl = sourceFileBS fl .| parseFastqC

sinkFastq :: (MonadResource m, MonadThrow m)
          => FilePath -> ConduitT Fastq o m ()
sinkFastq fl = mapC fastqToByteString .| unlinesAsciiC .| sinkFileBS fl

sinkFastqGzip :: (PrimMonad m, MonadThrow m, MonadResource m)
              => FilePath -> ConduitT Fastq o m ()
sinkFastqGzip fl = mapC fastqToByteString .| unlinesAsciiC .| gzip .| sinkFileBS fl

parseFastqC :: MonadThrow m => ConduitT B.ByteString Fastq m ()
parseFastqC = conduitParser fastqParser .| mapC snd
{-# INLINE parseFastqC #-}

fastqParser :: Parser Fastq
fastqParser = do
    _ <- skipWhile (/='@') >> char8 '@'
    ident <- A.takeTill isEndOfLine
    endOfLine
    sequence <- BS.filter (not . isEndOfLine) <$> takeTill (=='+')
    char8 '+' >> A.skipWhile (not . isEndOfLine) >> endOfLine
    score <- BS.filter (not . isEndOfLine) <$>
        A.scan 0 (f (B.length sequence))
    skipWhile (/='@')
    return $ Fastq ident sequence score
  where
    f n i x | i >= n = Nothing
            | isEndOfLine x = Just i
            | otherwise = Just $ i + 1
{-# INLINE fastqParser #-}

fastqToByteString :: Fastq -> B.ByteString
fastqToByteString (Fastq a b c) = "@" <> a <> "\n" <> b <> "\n+\n" <> c
{-# INLINE fastqToByteString #-}

-- | Get the mean and variance of quality scores at every position.
qualitySummary :: Monad m => ConduitT Fastq o m [(Double, Double)]
qualitySummary = mapC (map fromIntegral . decodeQualSc) .| meanVarianceC

meanVarianceC :: Monad m => ConduitT [Double] o m [(Double, Double)]
meanVarianceC = peekC >>= \case
    Nothing -> error "Empty input"
    Just x -> fst <$> foldlC f (replicate (length x) (0,0), 0)
  where
    f (acc, n) xs = let acc' = zipWith g acc xs in (acc', n')
      where
        n' = n + 1
        g (m, s) x = (m', s')
          where
            m' = m + d / fromIntegral n'
            s' = s + d * (x - m')
            d  = x - m
{-# INLINE meanVarianceC #-}

decodeQualSc :: Fastq -> [Int]
decodeQualSc = map (fromIntegral . (\x -> x - 33)) . BS.unpack .fastqSeqQual
{-# INLINE decodeQualSc #-}

pError :: Int -> Double
pError x = 10 ** (negate (fromIntegral x) / 10)
{-# INLINE pError #-}

{-
mkFastqRecord l1 l2 l3 l4 = Fastq (parseLine1 l1) (parseLine2 l2)
    (parseLine3 l3) (parseLine4 l4)
  where
    parseLine1 x | B.head x == '@' = B.tail x
                 | otherwise = error $ "Parse line 1 failed: " ++ B.unpack x
    parseLine2 x | B.all f x = x
                 | otherwise = error $ "Parse line 2 failed: " ++ B.unpack x
      where
        f 'C' = True
        f 'G' = True
        f 'T' = True
        f 'A' = True
        f 'N' = True
        f _   = False
    parseLine3 x | B.head x == '+' = B.tail x
                 | otherwise = error $ "Parse line 3 failed: " ++ B.unpack x
    parseLine4 x | BS.all f x = x
                 | otherwise = error $ "Parse line 4 failed: " ++ B.unpack x
      where
        f b = let b' = fromIntegral b :: Int
              in b' >= 33 && b' <= 126
-}

-- | Remove trailing 'A'
trimPolyA :: Int -> Fastq -> Fastq
trimPolyA n f@(Fastq a b c)
    | B.length trailing >= n = Fastq a b' $ B.take (B.length b') c
    | otherwise = f
  where
    (b', trailing) = B.spanEnd (=='A') b
