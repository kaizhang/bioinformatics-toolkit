{-# LANGUAGE BangPatterns          #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables   #-}
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
    , nucleotideFreq
    ) where

import qualified Data.ByteString.Char8 as B
import           Data.Char8            (toUpper)
import qualified Data.HashMap.Strict   as M
import qualified Data.HashSet          as S
import           Data.Proxy            (Proxy (..))
import           Prelude               hiding (length)

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

instance Semigroup (DNA a) where
    (<>) (DNA x) (DNA y) = DNA (x <> y)

instance Monoid (DNA a) where
    mempty = DNA B.empty
    mconcat dnas = DNA . B.concat . map toBS $ dnas

class BioSeq' s where
    toBS :: s a -> B.ByteString
    unsafeFromBS :: B.ByteString -> s a

    slice :: Int -> Int -> s a -> s a

    length :: s a -> Int
    length = B.length . toBS
    {-# MINIMAL toBS, slice, unsafeFromBS #-}

instance BioSeq' DNA where
    toBS (DNA s) = s
    unsafeFromBS = DNA
    slice i l (DNA s) = DNA . B.take l . B.drop i $ s

instance BioSeq' RNA where
    toBS (RNA s) = s
    unsafeFromBS = RNA
    slice i l (RNA s) = RNA . B.take l . B.drop i $ s

instance BioSeq' Peptide where
    toBS (Peptide s) = s
    unsafeFromBS = Peptide
    slice i l (Peptide s) = Peptide . B.take l . B.drop i $ s

class BioSeq' seq => BioSeq seq alphabet where
    alphabet :: Proxy (seq alphabet) -> S.HashSet Char
    fromBS :: B.ByteString -> Either String (seq alphabet)
    fromBS input = case B.mapAccumL fun Nothing input of
        (Nothing, r) -> Right $ unsafeFromBS r
        (Just e, _)  -> Left $ "Bio.Seq.fromBS: unknown character: " ++ [e]
      where
        fun (Just e) x = (Just e, x)
        fun Nothing x = let x' = toUpper x
                        in if x' `S.member` alphabet (Proxy :: Proxy (seq alphabet))
                            then (Nothing, x')
                            else (Just x', x')
    {-# MINIMAL alphabet #-}

instance BioSeq DNA Basic where
    alphabet _ = S.fromList "ACGT"

instance BioSeq DNA IUPAC where
    alphabet _ = S.fromList "ACGTNVHDBMKWSYR"

instance BioSeq DNA Ext where
    alphabet _ = undefined

instance BioSeq RNA Basic where
    alphabet _ = S.fromList "ACGU"

-- | O(n) Reverse complementary of DNA sequence.
rc :: DNA alphabet -> DNA alphabet
rc (DNA s) = DNA . B.map f . B.reverse $ s
  where
    f x = case x of
        'A' -> 'T'
        'C' -> 'G'
        'G' -> 'C'
        'T' -> 'A'
        _   -> x

-- | O(n) Compute GC content.
gcContent :: DNA alphabet -> Double
gcContent = (\(a,b) -> a / fromIntegral b) . B.foldl' f (0.0,0::Int) . toBS
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
                _   -> x + 0.5     -- "NMKYR"
        in (x', n+1)

-- | O(n) Compute single nucleotide frequency.
nucleotideFreq :: forall a . BioSeq DNA a => DNA a -> M.HashMap Char Int
nucleotideFreq dna = B.foldl' f m0 . toBS $ dna
  where
    m0 = M.fromList . zip (S.toList $ alphabet (Proxy :: Proxy (DNA a))) . repeat $ 0
    f m x = M.adjust (+1) x m
{-# INLINE nucleotideFreq #-}
