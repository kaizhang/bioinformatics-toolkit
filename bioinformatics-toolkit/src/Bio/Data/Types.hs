{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE OverloadedStrings #-}

module Bio.Data.Types where
    ( FileFormat(..)
    , File(..)
    , emptyFile
    , location
    , replication
    , format
    , keywords
    , Experiment(..)
    , eid
    , control
    , celltype
    , target
    , files
    , info
    ) where

import Data.Aeson
import Data.Aeson.TH (deriveJSON, defaultOptions)
import Shelly hiding (FilePath)
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as BC
import qualified Data.Text as T
import Control.Lens (makeFields, (^.), (.~))
import qualified Data.HashMap.Strict as M
import Crypto.Hash.MD5 (hash)
import Numeric (showHex)

data FileFormat = BamFile
                | BaiFile
                | BedFile
                | BedGZip
                | FastqFile
                | FastqGZip
                | BedgraphFile
                | BigWigFile
                | NarrowPeakFile
                | BroadPeakFile
                | Other
    deriving (Show, Read, Eq)

data File = File
    { fileLocation :: !FilePath
    , fileReplication :: !Int
    , fileFormat :: !FileFormat
    , fileKeywords :: ![String]
    , fileInfo :: !(M.HashMap String String)
    } deriving (Show, Read, Eq)

makeFields ''File

emptyFile :: File
emptyFile = File
    { fileLocation = ""
    , fileReplication = 0
    , fileFormat = Other
    , fileKeywords = []
    , fileInfo = M.empty
    }

data Experiment = Experiment
    { experimentEid :: !String
    , experimentCelltype :: !String
    , experimentTarget :: !String
    , experimentFiles :: ![File]
    , experimentInfo :: !(M.HashMap String String)
    , experimentControl :: !(Maybe String)
    } deriving (Show, Read, Eq)

makeFields ''Experiment

deriveJSON defaultOptions ''Format

instance FromJSON File where
    parseJSON (Object obj) = do
        path <- obj .: "path"
        File <$> return path <*>
                 obj .:? "rep" .!= 1 <*>
                 obj .:? "format" .!= getFormat path <*>
                 return [] <*>
                 return M.empty

instance FromJSON Experiment where
    parseJSON (Object obj) = do
        fls <- obj .: "files"
        let eid' = concat . map (flip showHex "") $ B.unpack $ hash $ BC.pack $
                   unlines $ map (^.location) fls
        Experiment <$> obj .:? "id" .!= eid' <*>
                       obj .:? "celltype" .!= "" <*>
                       obj .: "target" <*>
                       return fls <*>
                       return M.empty <*>
                       obj .:? "control"

getFormat :: FilePath -> Format
getFormat fl = case suf of
    "bam" -> Bam
    "bed" -> Bed
    "gz" -> case (snd $ T.breakOnEnd "." $ T.init pre) of
        "bed" -> BedGZip
        _ -> error "Unknown file format"
    _ -> error "Unknown file format"
  where
    (pre, suf) = T.breakOnEnd "." $ T.toLower $ T.pack fl
