module Bio.Utils.Misc (
      readInt
    , readDouble
    , bins
    , binBySize
) where

import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.ByteString.Lex.Double as L
import qualified Data.ByteString.Lex.Lazy.Double as LL
import Data.Maybe

class ToNum a where
    readInt :: a -> Int
    readDouble :: a -> Double

instance ToNum B.ByteString where
    readInt x = fst . fromMaybe raiseError . B.readInt $ x
      where raiseError = error $ "Fail to cast ByteString to Int:" ++ show x
    {-# INLINE readInt #-}

    readDouble x = fst . fromMaybe raiseError. L.readDouble $ x
      where raiseError = error $ "Fail to cast ByteString to Double:" ++ show x
    {-# INLINE readDouble #-}

instance ToNum BL.ByteString where
    readInt x = fst. fromMaybe raiseError. BL.readInt $ x
      where raiseError = error $ "Fail to cast ByteString to Int:" ++ show x
    {-# INLINE readInt #-}

    readDouble x = fst . fromMaybe raiseError . LL.readDouble $ x
      where raiseError = error $ "Fail to cast ByteString to Double:" ++ show x
    {-# INLINE readDouble #-}

-- | divide a given half-close-half-open region into fixed size 
-- half-close-half-open intervals, discarding leftovers
binBySize :: Int -> (Int, Int) -> [(Int, Int)]
binBySize step (start, end) = let xs = [start, start + step .. end]
                              in zip xs . tail $ xs
{-# INLINE binBySize #-}

-- | divide a given region into k equal size sub-regions, discarding leftovers
bins :: Int -> (Int, Int) -> [(Int, Int)]
bins binNum (start, end) = let k = (end - start) `div` binNum
                           in take binNum . binBySize k $ (start, end)
{-# INLINE bins #-}
