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
import Control.Monad.State.Lazy
import qualified Data.ByteString.Char8 as B
import Data.Conduit
import qualified Data.Conduit.List as CL
import System.IO

class FastaLike f where
    fromFastaRecord :: ( B.ByteString    -- ^ record header
                       , [B.ByteString]  -- ^ record body
                       )
                    -> f

    readFasta :: FilePath -> Source IO f
    readFasta fl = mapOutput fromFastaRecord . fastaReader $ fl

    -- | non-stream version, read whole file in memory
    readFasta' :: FilePath -> IO [f]
    readFasta' fl = readFasta fl $$ CL.consume
--    writeFasta :: FilePath -> [f] -> IO ()
    {-# MINIMAL fromFastaRecord #-}

instance BioSeq s a => FastaLike (s a) where
    fromFastaRecord (_, xs) = fromBS . B.concat $ xs
    {-# INLINE fromFastaRecord #-}

instance FastaLike Motif where
    fromFastaRecord (name, mat) = Motif name (toPWM mat)
    {-# INLINE fromFastaRecord #-}

fastaReader :: FilePath -> Source IO (B.ByteString, [B.ByteString])
fastaReader fl = do handle <- liftIO $ openFile fl ReadMode
                    loop [] handle
  where
    loop :: [B.ByteString] -> Handle -> Source IO (B.ByteString, [B.ByteString])
    loop acc h = do 
        eof <- liftIO $ hIsEOF h
        if eof then liftIO (hClose h) >> output acc else do
            l <- liftIO $ B.hGetLine h
            case () of
                _ | B.null l -> loop acc h  -- empty line, go to next line
                  | B.head l == '>' -> output acc >> loop [B.tail l] h
                  | otherwise -> loop (acc ++ [l]) h
    output (x:xs) = yield (x, xs)
    output _ = return ()
{-# INLINE fastaReader #-}
