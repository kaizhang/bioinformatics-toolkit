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
    , length
    -- * DNA related functions
    , rc
    -- * IO
    , Genome
    , GenomeH
    , gHOpen
    , gHClose
    , pack
    , getSeqs
    , getSeq
    , getIndex
    ) where

import Prelude hiding (length)
import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Lazy as M
import qualified Data.HashSet as S
import System.IO
import Data.List.Split
import Bio.Utils.Misc (readInt)
import Data.Char8
import Control.Monad
import Data.Monoid

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

instance Show (DNA a) where
    show (DNA s) = B.unpack s

instance Monoid (DNA a) where
    mempty = DNA B.empty
    mappend (DNA x) (DNA y) = DNA (x <> y)

class BioSeq' s where
    toBS :: s a -> B.ByteString
    length :: s a -> Int

instance BioSeq' DNA where
    toBS (DNA s) = s
    length = B.length . toBS

instance BioSeq' RNA where
    toBS (RNA s) = s
    length = B.length . toBS

instance BioSeq' Peptide where
    toBS (Peptide s) = s
    length = B.length . toBS

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

newtype Genome = G FilePath

newtype GenomeH = GH Handle

pack :: FilePath -> Genome
pack = G

gHOpen :: Genome -> IO GenomeH
gHOpen (G fl) = do h <- openFile fl ReadMode
                   return $ GH h

gHClose :: GenomeH -> IO ()
gHClose (GH h) = hClose h

type IndexTable = M.HashMap B.ByteString Int

-- | the number of characters before the start of genome
type OffSet = Int

type Query = (B.ByteString, Int, Int) -- (chr, start, end), zero based

getSeqs :: BioSeq s a => Genome -> [Query] -> IO [s a]
getSeqs g querys = do gH <- gHOpen g
                      (index, offset) <- getIndex gH
                      r <- mapM (getSeq gH index offset) querys
                      gHClose gH
                      return r

getSeq :: BioSeq s a => GenomeH -> IndexTable -> OffSet -> Query -> IO (s a)
getSeq (GH h) index offset (chr, start, end) = do 
    hSeek h AbsoluteSeek (fromIntegral pos)
    liftM fromBS $ B.hGet h (end - start + 1)
  where
    pos = offset + chrStart + start
    chrStart = M.lookupDefault (error $ "Bio.Seq.getSeq: Cannot find " ++ show chr) chr index

getIndex :: GenomeH -> IO (IndexTable, OffSet)
getIndex (GH h) = do header <- B.hGetLine h
                     return ( M.fromList . map f . chunksOf 2 . B.words $ header
                            , B.length header + 1 )
  where
    f [k, v] = (k, readInt v)
    f _ = error "error"
