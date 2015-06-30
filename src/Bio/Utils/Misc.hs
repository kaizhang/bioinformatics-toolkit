module Bio.Utils.Misc
    ( readInt
    , readDouble
    , bins
    , binBySize
    , binBySizeLeft
    ) where

import Data.ByteString.Char8 (ByteString)
import Data.ByteString.Lex.Fractional (readSigned, readExponential)
import Data.ByteString.Lex.Integral (readDecimal)
import Data.Maybe (fromMaybe)

readInt :: ByteString -> Int
readInt x = fst . fromMaybe errMsg . readSigned readDecimal $ x
  where
    errMsg = error $ "readInt: Fail to cast ByteString to Int:" ++ show x
{-# INLINE readInt #-}

readDouble :: ByteString -> Double
readDouble x = fst . fromMaybe errMsg . readSigned readExponential $ x
  where
    errMsg = error $ "readDouble: Fail to cast ByteString to Double:" ++ show x
{-# INLINE readDouble #-}

-- | divide a given half-close-half-open region into fixed size
-- half-close-half-open intervals, discarding leftovers
binBySize :: Int -> (Int, Int) -> [(Int, Int)]
binBySize step (start, end) = let xs = [start, start + step .. end]
                              in zip xs . tail $ xs
{-# INLINE binBySize #-}

-- | Including leftovers, the last bin may not have desired size.
binBySizeLeft :: Int -> (Int, Int) -> [(Int, Int)]
binBySizeLeft step (start, end) = let xs = [start, start + step .. end-1] ++ [end]
                                  in zip xs . tail $ xs
{-# INLINE binBySizeLeft #-}

-- | divide a given region into k equal size sub-regions, discarding leftovers
bins :: Int -> (Int, Int) -> [(Int, Int)]
bins binNum (start, end) = let k = (end - start) `div` binNum
                           in take binNum . binBySize k $ (start, end)
{-# INLINE bins #-}
