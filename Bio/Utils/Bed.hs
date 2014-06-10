{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DeriveGeneric #-}

module Bio.Utils.Bed (
      readInt
    , readDouble
    , bins
    , binBySize
    , BED(..)
    , chrom
    , chromStart
    , chromEnd
    , name
    , score
    , strand
    , readBED
    , writeBED 
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
import System.IO

-- | the type for BED format, as described in http://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
data BED = BED
    { _chrom :: B.ByteString
    , _chromStart :: {-# UNPACK #-} !Int
    , _chromEnd :: {-# UNPACK #-} !Int
    , _name :: Maybe B.ByteString
    , _score :: Maybe Double
    , _strand :: Maybe Bool  -- ^ True: "+", False: "-"
    } deriving (Read, Show, Generic)

makeLenses ''BED

instance Default BED

-- | FIXME: use handler
readBED :: FilePath -> IO [BED]
{-# INLINE readBED #-}
readBED fl = do content <- B.readFile fl
                return.map fromLine.B.lines $ content

writeBED :: FilePath -> [BED] -> IO ()
{-# INLINE writeBED #-}
writeBED fl beds = withFile fl WriteMode $ \h -> mapM_ (B.hPutStrLn h.toLine) beds

fromLine :: B.ByteString -> BED
{-# INLINE fromLine #-}
fromLine l = evalState (f (B.split '\t' l)) 1
  where
    f :: [B.ByteString] -> State Int BED
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

toLine :: BED -> B.ByteString
{-# INLINE toLine #-}
toLine (BED f1 f2 f3 f4 f5 f6) = B.intercalate "\t" [ f1
                                                    , (B.pack.show) f2
                                                    , (B.pack.show) f3
                                                    , fromMaybe "." f4
                                                    , score'
                                                    , strand'
                                                    ]
  where
    strand' | f6 == Just True = "+"
            | f6 == Just False = "-"
            | otherwise = "."
    score' = case f5 of
                 Just x -> (B.pack.show) x
                 _ -> "."

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
