{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE UndecidableInstances #-}
--------------------------------------------------------------------------------
-- |
-- Module      :  $Header$
-- Copyright   :  (c) 2014 Kai Zhang
-- License     :  MIT

-- Maintainer  :  kai@kzhang.org
-- Stability   :  experimental
-- Portability :  portable

-- Functions for processing Fasta files
--------------------------------------------------------------------------------

module Bio.Data.Fasta
    ( FastaLike(..)
    , fastaReader
    ) where

import Bio.Motif
import Bio.Seq
import qualified Data.ByteString.Char8 as B
import Conduit

class FastaLike f where
    fromFastaRecord :: ( B.ByteString    -- ^ record header
                       , [B.ByteString]  -- ^ record body
                       )
                    -> f

    readFasta :: FilePath -> Source (ResourceT IO) f
    readFasta fl = fastaReader fl =$= mapC fromFastaRecord

    -- | non-stream version, read whole file in memory
    readFasta' :: FilePath -> IO [f]
    readFasta' fl = runResourceT $ readFasta fl $$ sinkList
--    writeFasta :: FilePath -> [f] -> IO ()
    {-# MINIMAL fromFastaRecord #-}

instance BioSeq s a => FastaLike (s a) where
    fromFastaRecord (_, xs) = fromBS . B.concat $ xs
    {-# INLINE fromFastaRecord #-}

instance FastaLike Motif where
    fromFastaRecord (name, mat) = Motif name (toPWM mat)
    {-# INLINE fromFastaRecord #-}

fastaReader :: FilePath -> Source (ResourceT IO) (B.ByteString, [B.ByteString])
fastaReader fl = sourceFile fl =$= linesUnboundedAsciiC =$= loop []
  where
    loop acc = do
        x <- await
        case x of
            Just l -> case () of
                _ | B.null l -> loop acc  -- empty line, go to next line
                  | B.head l == '>' -> output (reverse acc) >> loop [B.tail l]
                  | otherwise -> loop (l:acc)
            Nothing -> output $ reverse acc
    output (x:xs) = yield (x, xs)
    output _ = return ()
{-# INLINE fastaReader #-}
