{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DeriveGeneric #-}

module Bio.RealWorld.Reactome
    ( getPathways
    ) where

import Data.Aeson
import GHC.Generics (Generic)
import qualified Data.Text as T
import           Network.HTTP.Simple

base :: String
base = "https://reactome.org/ContentService"

data Obj = Obj
    { className	:: Maybe T.Text
    , dbId :: Int 
    , displayName :: T.Text
    , schemaClass :: Maybe T.Text
    , stId :: Maybe T.Text
    , stIdVersion :: Maybe T.Text
    } deriving (Show, Generic)

instance ToJSON Obj
instance FromJSON Obj
 
-- | All Reactome top level pathways
getPathways :: String -> IO [Obj]
getPathways species = do
    req <- parseRequest url
    response <- httpJSON req
    return $ getResponseBody response
  where
    url = base ++ "/data/pathways/top/" ++ species

{-
pathwayAnalysis :: [B.ByteString]   -- ^ A list of identifiers
                -> IO B.ByteString
pathwayAnalysis ids = do
    initReq <- parseRequest base
    let request = urlEncodedBody [] initReq
    manager <- newManager tlsManagerSettings
    r <- fmap M.fromList $ runResourceT $ do
        response <- http request manager
        runConduit $ responseBody response .| linesUnboundedAsciiC .|
            (dropC 1 >> mapC ((\[a,b] -> (a,b)) . B.split '\t')) .| sinkList
    return $ map (flip M.lookup r) ids
  where
    params =
-}