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
    , BioSeq'(..)
    , BioSeq(..)
    -- * DNA related functions
    , rc
    ) where

import Prelude hiding (length)
import qualified Data.ByteString.Char8 as B
import qualified Data.HashSet as S
import Data.Char8
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
    mconcat dnas = DNA . B.concat . map toBS $ dnas

class BioSeq' s where
    toBS :: s a -> B.ByteString

    slice :: s a -> Int -> Int -> s a

    length :: s a -> Int
    length = B.length . toBS
    {-# MINIMAL toBS, slice #-}


instance BioSeq' DNA where
    toBS (DNA s) = s
    slice (DNA s) i l = DNA . B.take l . B.drop i $ s

instance BioSeq' RNA where
    toBS (RNA s) = s
    slice (RNA s) i l = RNA . B.take l . B.drop i $ s

instance BioSeq' Peptide where
    toBS (Peptide s) = s
    slice (Peptide s) i l = Peptide . B.take l . B.drop i $ s

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
