{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}
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
import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString as BS
import qualified Data.Attoparsec.ByteString as A
import Data.Word8
import Data.Attoparsec.ByteString.Char8
import Data.Conduit.Attoparsec
import GHC.Generics (Generic)
import Control.DeepSeq

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
    } deriving (Show, Eq, Generic, NFData)

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

data FQBuilder = Init FQBuilder
               | FQ1 B.ByteString FQBuilder
               | FQ2 B.ByteString FQBuilder
               | FQ3 B.ByteString FQBuilder
               | Complete

fqBuilder :: FQBuilder -> Fastq
fqBuilder = go ([], [], [])
  where
    go acc (Init bldr) = go acc bldr
    go (f1,f2,f3) (FQ1 x bldr) = go (x:f1, f2, f3) bldr
    go (f1,f2,f3) (FQ2 x bldr) = go (f1, x:f2, f3) bldr
    go (f1,f2,f3) (FQ3 x bldr) = go (f1, f2, x:f3) bldr
    go (f1,f2,f3) Complete = Fastq (B.concat $ reverse f1)
        (B.concat $ reverse f2) (B.concat $ reverse f3)
{-# INLINE fqBuilder #-}

parseFastqC :: Monad m => ConduitT B.ByteString Fastq m ()
parseFastqC = read1 >>= loop Init 'a'
  where
    read1 = await >>= \case
        Nothing -> error ""
        Just x -> return x
    loop acc st input = case st of
        'a' -> if B.head input == '@'
            then let (x, rest) = B.break (=='\n') $ B.tail input
                 in loop (FQ1 x) 'A' rest
            else error "Record does not start with \'@\'"
        'A' -> if B.null input
            then do
                input' <- read1
                let (x, rest) = B.break (=='\n') input'
                loop (acc . FQ1 x) 'A' rest
            else 
              let (x, rest) = B.break (=='\n') $ B.tail input in loop (acc . FQ2 x) 'b' rest
        'b' -> if B.null input
            then do
                input' <- read1
                let (x, rest) = B.break (=='\n') input'
                loop (acc . FQ2 x) 'b' rest
            else 
              if (B.head $ B.tail input) == '+'
                  then loop acc 'c' $ B.drop 2 input
                  else 
                      let (x, rest) = B.break (=='\n') $ B.tail input
                      in loop (acc . FQ2 x) 'b' rest
        'c' -> let (x, rest) = B.break (=='\n') input
               in if B.null rest
                  then do
                      bs <- read1
                      loop acc 'c' bs 
                  else 
                      loop acc 'd' $ B.tail rest
        'd' -> do
            input' <- if B.null input then read1 else return input
            let (x, rest) = B.break (=='\n') input'
            if B.null rest 
                then loop (acc . FQ3 x) 'd' rest
                else loop (acc . FQ3 x) 'D' $ B.tail rest
        'D' -> do 
            if B.null input
              then await >>= \case
                  Nothing -> yield $ fqBuilder $ acc Complete
                  Just input' -> if B.head input' == '@'
                      then do
                          yield $ fqBuilder $ acc Complete
                          let (x, rest) = B.break (=='\n') $ B.tail input'
                          loop (Init . FQ1 x) 'A' rest
                      else let (x, rest) = B.break (=='\n') input' in loop (acc . FQ3 x) 'D' $ B.tail rest
            else if B.head input == '@'
                then do
                    yield $ fqBuilder $ acc Complete
                    let (x, rest) = B.break (=='\n') $ B.tail input
                    loop (Init . FQ1 x) 'A' rest
                else let (x, rest) = B.break (=='\n') input in loop (acc . FQ3 x) 'D' $ B.tail rest
{-# INLINE parseFastqC #-}

{-
parseFastqC :: MonadThrow m => ConduitT B.ByteString Fastq m ()
parseFastqC = conduitParser fastqParser .| mapC snd
{-# INLINE parseFastqC #-}

fastqParser :: Parser Fastq
fastqParser = do
    _ <- skipWhile (/='@') >> char8 '@'
    ident <- A.takeTill isEndOfLine
    endOfLine
    sequ <- BS.filter (not . isEndOfLine) <$> takeTill (=='+')
    char8 '+' >> A.skipWhile (not . isEndOfLine) >> endOfLine
    score <- BS.filter (not . isEndOfLine) <$>
        A.scan 0 (f (B.length sequ))
    skipWhile (/='@')
    return $ Fastq ident sequ score
  where
    f n i x | i >= n = Nothing
            | isEndOfLine x = Just i
            | otherwise = Just $ i + 1
{-# INLINE fastqParser #-}
-}

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