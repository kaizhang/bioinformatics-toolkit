{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Bio.Seq 
    ( 
    -- * Alphabet
      Basic
    , IUPAC
    , Ext
    -- * Sequence types
    , DNA
    , RNA
    , Peptide
    , BioSeq (..)
    , toBS
    -- * DNA related functions
    , rc
    -- * IO
    , getSeqs
    , getSeq
    , getIndex
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Lazy as M
import qualified Data.HashSet as S
import System.IO
import Data.List.Split
import Bio.Utils.Bed (readInt)
import Data.Char8
import Control.Monad

-- | Alphabet defined by http://www.chem.qmul.ac.uk/iupac/
-- | Standard unambiguous alphabet
data Basic

-- | full IUPAC alphabet, including ambiguous letters
data IUPAC

-- | extended alphabet
data Ext

-- | DNA sequence
newtype DNA alphabet = DNA B.ByteString

-- | RNA sequence
newtype RNA alphabet = RNA B.ByteString

-- | Peptide sequence
newtype Peptide alphabet = Peptide B.ByteString

class BioSeq' s where
    toBS :: s a -> B.ByteString

instance BioSeq' DNA where
    toBS (DNA s) = s

instance BioSeq' RNA where
    toBS (RNA s) = s

instance BioSeq' Peptide where
    toBS (Peptide s) = s

class BioSeq' s => BioSeq s a where
    fromBS :: B.ByteString -> s a 

instance BioSeq DNA Basic where
    fromBS = DNA . B.map (f.toUpper)
      where
        f x | x `S.member` alphabet = x
            | otherwise = error "error"
        alphabet = S.fromList "ACGT"

instance BioSeq DNA IUPAC where
    fromBS = DNA . B.map (f.toUpper)
      where
        f x | x `S.member` alphabet = x
            | otherwise = error "error"
        alphabet = S.fromList "ACGTNVHDBMKWSYR"

instance BioSeq DNA Ext where
    fromBS = undefined

instance BioSeq RNA Basic where
    fromBS = RNA . B.map (f.toUpper)
      where
        f x | x `S.member` alphabet = x
            | otherwise = error "error"
        alphabet = S.fromList "ACGU"

-- | reverse complementary of DNA sequence
rc :: DNA alphabet -> DNA alphabet
rc (DNA s) = DNA . B.map f . B.reverse $ s
  where
    f x = case x of
        'A' -> 'T'
        'C' -> 'G'
        'G' -> 'C'
        'T' -> 'A'
        _ -> x

type IndexTable = M.HashMap B.ByteString Int

type OffSet = Int

type Query = (B.ByteString, Int, Int) -- (chr, start, end), zero based

getSeqs :: BioSeq s a => FilePath -> [Query] -> IO [s a]
getSeqs fl querys = do h <- openFile fl ReadMode
                       (index, offset) <- getIndex h
                       r <- mapM (getSeq h index offset) querys
                       hClose h
                       return r

getSeq :: BioSeq s a => Handle -> IndexTable -> OffSet -> Query -> IO (s a)
getSeq h index offset (chr, start, end) = do 
    hSeek h AbsoluteSeek (fromIntegral pos)
    liftM fromBS $ B.hGet h (end - start + 1)
  where
    pos = offset + chrStart + start
    chrStart = M.lookupDefault (error $ "Bio.Seq.getSeq: Cannot find " ++ show chr) chr index

getIndex :: Handle -> IO (IndexTable, OffSet)
getIndex h = do header <- B.hGetLine h
                return (M.fromList.map f.chunksOf 2.B.words $ header, B.length header + 1)
  where
    f [k, v] = (k, readInt v)
    f _ = error "error"
