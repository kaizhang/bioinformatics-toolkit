{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE BangPatterns #-}

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
    , gcContent
    ) where

import Prelude hiding (length)
import qualified Data.ByteString.Char8 as B
import qualified Data.HashSet as S
import Data.Char8 (toUpper)
import Data.Monoid (Monoid(..))

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
    mappend (DNA x) (DNA y) = DNA (x `mappend` y)
    mconcat dnas = DNA . B.concat . map toBS $ dnas

class BioSeq' s where
    toBS :: s a -> B.ByteString

    slice :: Int -> Int -> s a -> s a

    length :: s a -> Int
    length = B.length . toBS
    {-# MINIMAL toBS, slice #-}


instance BioSeq' DNA where
    toBS (DNA s) = s
    slice i l (DNA s) = DNA . B.take l . B.drop i $ s

instance BioSeq' RNA where
    toBS (RNA s) = s
    slice i l (RNA s) = RNA . B.take l . B.drop i $ s

instance BioSeq' Peptide where
    toBS (Peptide s) = s
    slice i l (Peptide s) = Peptide . B.take l . B.drop i $ s

class BioSeq' s => BioSeq s a where
    alphabet :: s a -> S.HashSet Char
    fromBS :: B.ByteString -> s a 

instance BioSeq DNA Basic where
    alphabet _ = S.fromList "ACGT"
    fromBS = DNA . B.map (f . toUpper)
      where
        f x | x `S.member` alphabet (undefined :: DNA Basic) = x
            | otherwise = error $ "Bio.Seq.fromBS: unknown character: " ++ [x]

instance BioSeq DNA IUPAC where
    alphabet _ = S.fromList "ACGTNVHDBMKWSYR"
    fromBS = DNA . B.map (f . toUpper)
      where
        f x | x `S.member` alphabet (undefined :: DNA IUPAC) = x
            | otherwise = error $ "Bio.Seq.fromBS: unknown character: " ++ [x]

instance BioSeq DNA Ext where
    alphabet = undefined
    fromBS = undefined

instance BioSeq RNA Basic where
    alphabet _ = S.fromList "ACGU"
    fromBS = RNA . B.map (f . toUpper)
      where
        f x | x `S.member` alphabet (undefined :: RNA Basic) = x
            | otherwise = error $ "Bio.Seq.fromBS: unknown character: " ++ [x]

-- | O(n). reverse complementary of DNA sequence
rc :: DNA alphabet -> DNA alphabet
rc (DNA s) = DNA . B.map f . B.reverse $ s
  where
    f x = case x of
        'A' -> 'T'
        'C' -> 'G'
        'G' -> 'C'
        'T' -> 'A'
        _ -> x

-- | O(n). compute GC content
gcContent :: DNA alphabet -> Double
gcContent = (\(a,b) -> a / fromIntegral b) . B.foldl f (0.0,0::Int) . toBS
  where
    f (!x,!n) c =
        let x' = case c of
                'A' -> x
                'C' -> x + 1
                'G' -> x + 1
                'T' -> x
                'H' -> x + 0.25
                'D' -> x + 0.25
                'V' -> x + 0.75
                'B' -> x + 0.75
                'S' -> x + 1
                'W' -> x
                _ -> x + 0.5     -- "NMKYR"
        in (x', n+1)
