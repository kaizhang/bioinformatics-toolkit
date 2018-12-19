-- NIH Genomic Data Commons
{-# LANGUAGE OverloadedStrings #-}

module Bio.RealWorld.GDC
    (downloadData) where

import           Conduit
import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B
import Data.Maybe (fromJust)
import           Network.HTTP.Conduit

baseurl :: String
baseurl = "https://api.gdc.cancer.gov/"

-- | Download data
downloadData :: String    -- ^ UUID
             -> FilePath  -- ^ Output dir
             -> IO FilePath
downloadData uuid dir = do
     request <- parseRequest url
     manager <- newManager tlsManagerSettings
     runResourceT $ do
         response <- http request manager
         let filename = T.unpack $ snd $ T.breakOnEnd "filename=" $ T.pack $
                B.unpack $ fromJust $ lookup "Content-Disposition" $
                responseHeaders response
         runConduit $ responseBody response .| sinkFileBS (dir ++ "/" ++ filename)
         return filename
  where
    url = baseurl ++ "data/" ++ uuid
{-# INLINE downloadData #-}