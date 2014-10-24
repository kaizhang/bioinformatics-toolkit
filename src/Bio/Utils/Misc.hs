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

-- | divide a given region into fixed size fragments
binBySize :: (Int, Int) -> Int -> [(Int, Int)]
binBySize (start, end) step =
    let binNum = ceiling $ fromIntegral (end - start + 1) / (fromIntegral step :: Double)
    in take binNum $ zip [start,start+step..] [start+step-1,start+2*step-1..]

-- | divide a given region into k equal size sub-regions
bins :: (Int, Int) -> Int -> [(Int, Int)]
bins (start, end) binNum = 
    let step = ceiling $ fromIntegral (end - start + 1) / (fromIntegral binNum :: Double)
    in take binNum $ zip [start,start+step..] [start+step-1,start+2*step-1..]
