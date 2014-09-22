{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE UndecidableInstances #-}

module Bio.Utils.Fasta
    ( FastaFormat(..)
    , fastaReader
    ) where

import Bio.Motif
import Bio.Seq
import Control.Monad.State.Lazy
import Data.Conduit
import qualified Data.ByteString.Char8 as B
import System.IO

class FastaFormat f where
    fromFastaRecord :: ( B.ByteString    -- ^ record header
                       , [B.ByteString]  -- ^ record body
                       )
                    -> f
    readFasta :: FilePath -> Source IO f
    readFasta fl = mapOutput fromFastaRecord . fastaReader $ fl
--    writeFasta :: FilePath -> [f] -> IO ()
    {-# MINIMAL fromFastaRecord #-}

instance BioSeq s a => FastaFormat (s a) where
    fromFastaRecord (_, xs) = fromBS . B.concat $ xs
    {-# INLINE fromFastaRecord #-}

instance FastaFormat Motif where
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
