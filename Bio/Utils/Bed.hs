{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE DeriveGeneric #-}

module Bio.Utils.Bed (
      BED(..)
    , chrom
    , chromStart
    , chromEnd
    , name
    , score
    , strand
    , fetchSeq
    , readBED
    , writeBED 
) where

import qualified Data.ByteString.Char8 as B
import Bio.Seq
import Bio.Utils.Misc (readInt, readDouble)
import Data.Maybe
import Data.Conduit
import Control.Lens
import Control.Monad.State.Strict
import Data.Default.Generics
import GHC.Generics
import System.IO

-- | the type for BED format, as described in http://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
data BED = BED
    { _chrom :: !B.ByteString
    , _chromStart :: {-# UNPACK #-} !Int
    , _chromEnd :: {-# UNPACK #-} !Int
    , _name :: !(Maybe B.ByteString)
    , _score :: !(Maybe Double)
    , _strand :: !(Maybe Bool)  -- ^ True: "+", False: "-"
    } deriving (Read, Show, Generic)

makeLenses ''BED

instance Default BED

readBED :: FilePath -> Source IO BED
readBED fl = do handle <- liftIO $ openFile fl ReadMode
                loop handle
  where
    loop h = do eof <- liftIO $ hIsEOF h
                if eof 
                   then liftIO $ hClose h
                   else do
                       line <- liftIO $ B.hGetLine h
                       yield $ fromLine line
                       loop h
{-# INLINE readBED #-}

fetchSeq :: BioSeq DNA a => Genome -> Conduit BED IO (DNA a)
fetchSeq g = do gH <- liftIO $ gHOpen g
                (table, offset) <- liftIO $ getIndex gH
                conduitWith gH table offset
                liftIO $ gHClose gH
  where
    conduitWith h index' offset' = do 
        bed <- await
        case bed of
            Just (BED chr start end _ _ isForward) -> do 
                dna <- liftIO $ getSeq h index' offset' (chr, start, end)
                case isForward of
                    Just False -> yield $ rc dna
                    _ -> yield dna
                conduitWith h index' offset'
            _ -> return ()
{-# INLINE fetchSeq #-}

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
