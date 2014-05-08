{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

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

-- | the type for BED format, as described in http://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
data BED = BED
    { _chrom :: BL.ByteString
    , _chromStart :: {-# UNPACK #-} !Int
    , _chromEnd :: {-# UNPACK #-} !Int
    , _name :: Maybe BL.ByteString
    , _score :: Maybe Double
    , _strand :: Maybe Bool
    } deriving (Show)

makeLenses ''BED

readBED :: FilePath -> IO [BED]
readBED fl = do content <- BL.readFile fl
                return.map fromLine.BL.lines $ content

fromLine :: BL.ByteString -> BED
{-# INLINE fromLine #-}
fromLine l = case BL.words l of
    [chr,s,e] -> BED chr (readInt s) (readInt e) Nothing Nothing Nothing
    [chr,s,e,nm,sc] -> BED chr (readInt s) (readInt e) (getName nm) (getScore sc) Nothing
    [chr,s,e,nm,sc,str] -> BED chr (readInt s) (readInt e) (getName nm) (getScore sc) (getStrand str)
    _ -> error "Read BED fail: Incorrect number of fields"
  where
    getName nm | nm == "." = Nothing
               | otherwise = Just nm
    getScore s | s == "." = Nothing
               | otherwise = Just . readDouble $ s
    getStrand str | str == "-" = Just False
                  | str == "+" = Just True
                  | otherwise = Nothing

class ToNum a where
    readInt :: a -> Int
    readDouble :: a -> Double

instance ToNum B.ByteString where
    readInt = fst.fromJust.B.readInt
    readDouble = fst . fromJust . L.readDouble

instance ToNum BL.ByteString where
    readInt = fst.fromJust.BL.readInt
    readDouble = fst . fromJust . LL.readDouble

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
