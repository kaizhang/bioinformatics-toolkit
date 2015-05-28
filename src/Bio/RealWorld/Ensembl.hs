{-# LANGUAGE OverloadedStrings #-}
module Bio.RealWorld.Ensembl
    ( lookup
    ) where

import Prelude hiding (lookup)
import Data.Aeson
import Data.List.Split (chunksOf)
import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict as M
import Network.HTTP.Conduit

import Bio.RealWorld.ID (BioID(..), EnsemblID)

base :: String
base = "http://rest.ensembl.org/"

lookup :: [EnsemblID] -> IO (Either String Object)
lookup xs = do
    rs <- mapM lookupHelp $ chunksOf 1000 xs
    return $ foldl1 f rs
  where
    f a b = do
        a' <- a
        b' <- b
        return $ M.union a' b'

lookupHelp :: [EnsemblID] -> IO (Either String Object)
lookupHelp xs = do
    initReq <- parseUrl url
    let request = initReq { method = "POST" 
                          , requestHeaders = [("Content-type", "application/json")]
                          , requestBody = body
                          }
    r <- withManager $ \manager -> httpLbs request manager
    return . eitherDecode . responseBody $ r
  where
    url = base ++ "/lookup/id/"
    ids = B.pack $ show $ map fromID xs
    body = RequestBodyBS $ B.intercalate "" ["{ \"ids\" :", ids, "}"]
{-# INLINE lookupHelp #-}
