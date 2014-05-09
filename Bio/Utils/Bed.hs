{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DeriveGeneric #-}

module Bio.Utils.Bed (
      readInt
    , readDouble
    , bins
    , binBySize
    , BED
    , chrom
    , chromStart
    , chromEnd
    , name
    , score
    , strand
    , readBED
) where

import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.ByteString.Lex.Double as L
import qualified Data.ByteString.Lex.Lazy.Double as LL
import Data.Maybe
import Control.Lens
import Control.Monad.State.Strict
import Data.Default.Generics
import GHC.Generics

-- | the type for BED format, as described in http://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
data BED = BED
    { _chrom :: BL.ByteString
    , _chromStart :: {-# UNPACK #-} !Int
    , _chromEnd :: {-# UNPACK #-} !Int
    , _name :: Maybe BL.ByteString
    , _score :: Maybe Double
    , _strand :: Maybe Bool
    } deriving (Show, Generic)

makeLenses ''BED

instance Default BED

readBED :: FilePath -> IO [BED]
readBED fl = do content <- BL.readFile fl
                return.map fromLine.BL.lines $ content

fromLine :: BL.ByteString -> BED
{-# INLINE fromLine #-}
fromLine l = evalState (f (BL.split '\t' l)) 1
  where
    f :: [BL.ByteString] -> State Int BED
    f [] = do i <- get
              if i <= 3 then error "Read BED fail: Incorrect number of fields"
                        else return def
    f (x:xs) = do 
        i <- get
        put (i+1)
        bed <- f xs
        case i of
            1 -> return $ chrom .~ x $ bed
            2 -> return $ chromStart .~ readInt x $ bed
            3 -> return $ chromEnd .~ readInt x $ bed
            4 -> return $ name .~ guard' x $ bed
            5 -> return $ score .~ getScore x $ bed
            6 -> return $ strand .~ getStrand x $ bed
            _ -> return def

    guard' x | x == "." = Nothing
             | otherwise = Just x
    getScore x | x == "." = Nothing
               | otherwise = Just . readDouble $ x
    getStrand str | str == "-" = Just False
                  | str == "+" = Just True
                  | otherwise = Nothing

class ToNum a where
    readInt :: a -> Int
    readDouble :: a -> Double

instance ToNum B.ByteString where
    readInt = fst.fromJust.B.readInt
    {-# INLINE readInt #-}
    readDouble = fst . fromJust . L.readDouble
    {-# INLINE readDouble #-}

instance ToNum BL.ByteString where
    readInt = fst.fromJust.BL.readInt
    {-# INLINE readInt #-}
    readDouble = fst . fromJust . LL.readDouble
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
