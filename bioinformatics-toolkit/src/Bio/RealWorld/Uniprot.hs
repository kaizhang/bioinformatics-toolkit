{-# LANGUAGE OverloadedStrings #-}

module Bio.RealWorld.Uniprot
    ( mapID
    ) where

import           Conduit
import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict   as M
import           Network.HTTP.Conduit

import           Bio.RealWorld.ID

base :: String
base = "http://www.uniprot.org/uploadlists/"

mapID :: [B.ByteString]   -- ^ A list of IDs
      -> B.ByteString     -- ^ From database
      -> B.ByteString     -- ^ To database
      -> IO [Maybe B.ByteString]
mapID ids from to = do
    initReq <- parseRequest base
    let request = setQueryString query initReq
            { method = "GET"
            , requestHeaders = [("User-Agent", "kk@test.org")]
            }
    manager <- newManager tlsManagerSettings
    r <- fmap M.fromList $ runResourceT $ do
        response <- http request manager
        runConduit $ responseBody response .| linesUnboundedAsciiC .|
            (dropC 1 >> mapC ((\[a,b] -> (a,b)) . B.split '\t')) .| sinkList
    return $ map (flip M.lookup r) ids
  where
    query = [ ("from", Just from)
            , ("to", Just to)
            , ("format", Just "tab")
            , ("query", Just $ B.unwords ids)
            ]
