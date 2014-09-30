{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DeriveGeneric #-}
--------------------------------------------------------------------------------
-- |
-- Module      :  $Header$
-- Description :  Data types for ENCODE record
-- Copyright   :  (c) Kai Zhang
-- License     :  MIT

-- Maintainer  :  kai@kzhang.org
-- Stability   :  experimental
-- Portability :  portable

-- Data types for ENCODE record
--------------------------------------------------------------------------------

module Bio.Data.ENCODE.Types
    ( Record(..)
    , parseRecord
    ) where

import Data.Aeson
import Data.Aeson.Types
import qualified Data.Text as T
import GHC.Generics

-- | represent the search results
data Record = Record
    { _dbxrefs :: ![T.Text]
    , _status :: !T.Text
    , _accession :: !T.Text
    , _hub :: !T.Text
    , _lab :: !T.Text
    , _systemSlims :: ![T.Text]
    , _monthReleased :: !T.Text
    , _organSlims :: ![T.Text]
    , _replicates :: ![T.Text]
    , _biosampleType :: !T.Text
    , _datasetType :: !T.Text
    , _uuid :: !T.Text
    , _references :: ![T.Text]
    , _assayTermName :: !T.Text
    , _aliases :: ![T.Text]
    , _originalFiles :: ![T.Text]
    , _runType :: !T.Text
    , _revokedFiles :: ![T.Text]
    , _developmentalSlims :: ![T.Text]
    , _schemaVersion :: !T.Text
    , _documents :: ![T.Text]
    , _synonyms :: ![T.Text]
    , _relatedFiles :: ![T.Text]
    , _files :: ![T.Text]
    , _submittedBy :: !T.Text
    , _id :: !T.Text
    , _type :: ![T.Text]
    , _alternateAccessions :: ![T.Text]
    , _biosampleTermName :: !T.Text
    , _dateCreated :: !T.Text
    , _dateReleased :: !T.Text
    , _assembly :: !T.Text
    , _possibleControls :: ![T.Text]
    , _description :: !T.Text
    , _target :: !T.Text
    , _award :: !T.Text
    , _biosampleTermId :: !T.Text
    } deriving (Show, Generic)

instance FromJSON Record

instance ToJSON Record

parseRecord :: Value -> Parser Record
parseRecord = genericParseJSON defaultOptions
                  { fieldLabelModifier = modify . tail}
  where
    modify x | x `elem` ["id", "type"] = camelTo '_' ('@':x)
             | otherwise = camelTo '_' x
