module Bio.Utils.BitVector
    ( BitVector
    , BitMVector
    , size
    , (!)
    , set
    , clear
    , unsafeFreeze
    , zeros
    , toList
    ) where

import qualified Data.Vector.Unboxed as U
import Control.Monad.Primitive
import qualified Data.Vector.Unboxed.Mutable as UM
import Data.Word
import Data.Bits
import Text.Printf (printf)

data BitVector = BitVector Int (U.Vector Word8)

data BitMVector s = BitMVector Int (UM.MVector s Word8)

size :: BitVector -> Int
size (BitVector n _) = n

(!) :: BitVector -> Int -> Bool
(!) = index

index :: BitVector -> Int -> Bool
index (BitVector n v) idx
    | idx >= n = error $ printf "index out of bounds (%d,%d)" idx n
    | otherwise = testBit (v `U.unsafeIndex` i) j
  where
    i = idx `div` 8
    j = idx `mod` 8

set :: PrimMonad m => BitMVector (PrimState m) -> Int -> m ()
set (BitMVector _ mv) idx = UM.modify mv ((flip setBit) j) i
  where
    i = idx `div` 8
    j = idx `mod` 8

clear :: PrimMonad m => BitMVector (PrimState m) -> Int -> m ()
clear (BitMVector _ mv) idx = UM.modify mv ((flip clearBit) j) i
  where
    i = idx `div` 8
    j = idx `mod` 8

unsafeFreeze :: PrimMonad m => BitMVector (PrimState m) -> m BitVector
unsafeFreeze (BitMVector n mv) = U.unsafeFreeze mv >>= return . BitVector n

zeros :: PrimMonad m => Int -> m (BitMVector (PrimState m))
zeros n = UM.replicate n' 0 >>= return . BitMVector n
  where
    n' = if j == 0 then i else i + 1
    i = n `div` 8
    j = n `mod` 8

toList :: BitVector -> [Bool]
toList bv = flip map [0..n-1] $ \i -> bv ! i
  where
    n = size bv