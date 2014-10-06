{-# LANGUAGE OverloadedStrings #-}
--------------------------------------------------------------------------------
-- |
-- Module      :  $Header$
-- Description :  Search and download data from ENCODE project
-- Copyright   :  (c) Kai Zhang
-- License     :  MIT

-- Maintainer  :  kai@kzhang.org
-- Stability   :  experimental
-- Portability :  portable

-- Search and download data from ENCODE project
--------------------------------------------------------------------------------

module Bio.Data.ENCODE
    ( search
    ) where

import Bio.Data.ENCODE.Types
import Data.Aeson
import Data.Aeson.Types
import Data.List
import Network.HTTP.Conduit

base :: String
base = "https://www.encodeproject.org/"

search :: [String] -> [String] -> IO [Either String Record]
search terms constraints = do 
    initReq <- parseUrl url
    let request = initReq { method = "GET" 
                          , requestHeaders = [("accept", "application/json")]
                          }
    r <- withManager $ \manager -> httpLbs request manager
    case eitherDecode (responseBody r) >>= parseEither parser of
        Left msg -> error msg
        Right x -> return x
  where
    url = intercalate "" [ base                   -- base url
                         , "search/?"
                         , mkSearchTerm terms
                         , "&"
                         , mkConstraint constraints
                         , "&frame=object"        -- get all object properties
                         ]

    mkSearchTerm [] = ""
    mkSearchTerm xs = "searchTerm=" ++ intercalate "+" xs

    mkConstraint [] = ""
    mkConstraint xs = intercalate "&" xs

    parser x = do xs <- withObject "ENCODE_JSON" (.: "@graph") x
                  return $ map (parseEither parseRecord) xs
