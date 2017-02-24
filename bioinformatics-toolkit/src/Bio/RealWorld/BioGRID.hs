{-# LANGUAGE OverloadedStrings #-}

module Bio.RealWorld.BioGRID
    ( TAB2(..)
    , fetchByGeneNames
    ) where

import Network.HTTP.Conduit
import Data.List
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T

accessKey :: String
accessKey = "accessKey=6168b8d02b2aa2e9a45af6f3afac4461"

base :: String
base = "http://webservice.thebiogrid.org/"

-- | BioGRID tab2 format
data TAB2 = TAB2
    { _biogridId :: B.ByteString
    , _entrezIdA :: B.ByteString
    , _entrezIdB :: B.ByteString
    , _biogridIdA :: B.ByteString
    , _biogridIdB :: B.ByteString
    , _systematicNameA :: T.Text
    , _systematicNameB :: T.Text
    , _symbolA :: T.Text
    , _symbolB :: T.Text
    , _synonymsA :: [T.Text]
    , _synonymsB :: [T.Text]
    , _experimentalSystemName :: T.Text
    , _experimentalSystemType :: T.Text
    , _firstAuthor :: T.Text
    , _pubmedId :: B.ByteString
    , _organismIdA :: B.ByteString
    , _organismIdB :: B.ByteString
    , _throughput :: T.Text
    , _score :: Maybe Double
    , _ptm :: T.Text
    , _phenotypes :: [T.Text]
    , _qualifications :: [T.Text]
    , _tags :: [T.Text]
    , _source :: T.Text
    } deriving (Show)

parseAsTab2 :: BL.ByteString -> TAB2
parseAsTab2 l = TAB2 (BL.toStrict $ xs!!0)
                     (BL.toStrict $ xs!!1)
                     (BL.toStrict $ xs!!2)
                     (BL.toStrict $ xs!!3)
                     (BL.toStrict $ xs!!4)
                     (T.pack $ BL.unpack $ xs!!5)
                     (T.pack $ BL.unpack $ xs!!6)
                     (T.pack $ BL.unpack $ xs!!7)
                     (T.pack $ BL.unpack $ xs!!8)
                     (T.splitOn "|" $ T.pack $ BL.unpack $ xs!!9)
                     (T.splitOn "|" $ T.pack $ BL.unpack $ xs!!10)
                     (T.pack $ BL.unpack $ xs!!11)
                     (T.pack $ BL.unpack $ xs!!12)
                     (T.pack $ BL.unpack $ xs!!13)
                     (BL.toStrict $ xs!!14)
                     (BL.toStrict $ xs!!15)
                     (BL.toStrict $ xs!!16)
                     (T.pack $ BL.unpack $ xs!!17)
                     (getScore $ BL.unpack $ xs!!18)
                     (T.pack $ BL.unpack $ xs!!19)
                     (T.splitOn "|" $ T.pack $ BL.unpack $ xs!!20)
                     (T.splitOn "|" $ T.pack $ BL.unpack $ xs!!21)
                     (T.splitOn "|" $ T.pack $ BL.unpack $ xs!!22)
                     (T.pack $ BL.unpack $ xs!!23)
  where
    xs = BL.split '\t' l
    getScore "-" = Nothing
    getScore x = Just $ read x

-- | retreive first 10,000 records
fetchByGeneNames :: [String] -> IO [TAB2]
fetchByGeneNames genes = do
    initReq <- parseRequest $ intercalate "&" [url, geneList, tax, accessKey]
    let request = initReq { method = "GET"
                          , requestHeaders = [("Content-type", "text/plain")]
                          }
    manager <- newManager tlsManagerSettings
    r <- httpLbs request manager
    return $ map parseAsTab2 $ BL.lines $ responseBody r
  where
    url = base ++ "/interactions/?searchNames=ture&includeInteractors=false"
    geneList = "geneList=" ++ intercalate "|" genes
    tax = "taxId=9606"
{-# INLINE fetchByGeneNames #-}
