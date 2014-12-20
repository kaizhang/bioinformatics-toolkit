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

module Bio.RealWorld.ENCODE
    ( search
    , showResult
    , (|@)
    , (|!)
    , as
    , (&)
    ) where

import Data.Aeson
import Data.Aeson.Types
import Data.Aeson.Encode.Pretty (encodePretty)
import qualified Data.HashMap.Lazy as M
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Sequence as S
import qualified Data.Foldable as F
import qualified Data.Text as T
import qualified Data.Vector as V
import Network.HTTP.Conduit
import Data.Default.Class
import Data.Monoid

data KeyWords = KeyWords (S.Seq String)  -- ^ terms
                         (S.Seq String)  -- ^ constraints

instance Default KeyWords where
    def = KeyWords S.empty $ S.fromList ["frame=object", "limit=all"]

instance Show KeyWords where
    show (KeyWords x y) = f x ++ g y
      where
        f x' | S.null x' = ""
             | otherwise = "searchTerm=" ++ F.foldr1 (\a b -> b ++ ('+':a)) x' ++ "&"
        g y' | S.null y' = ""
             | otherwise =  F.foldr1 (\a b -> b ++ ('&':a)) y'

instance Monoid KeyWords where
    mempty = KeyWords S.empty S.empty
    mappend (KeyWords a b) (KeyWords a' b') = KeyWords (a S.>< a') (b S.>< b')

base :: String
base = "https://www.encodeproject.org/"

-- | general search using keywords and a set of constraints. Example:
-- search ["chip", "sp1"] ["type=experiment"]
search :: KeyWords -> IO (Either String [Value])
search kw = do 
    initReq <- parseUrl url
    let request = initReq { method = "GET" 
                          , requestHeaders = [("accept", "application/json")]
                          }
    r <- withManager $ \manager -> httpLbs request manager
    return $ (eitherDecode . responseBody) r >>=
                parseEither (\x -> withObject "ENCODE_JSON" (.: "@graph") x)
  where
    url = base ++ "search/?" ++ show kw

showResult :: Value -> IO ()
showResult = B.putStrLn . encodePretty

-- * common search

terms :: [String] -> KeyWords
terms xs = KeyWords (S.fromList xs) S.empty

assayIs :: String -> KeyWords
assayIs x = KeyWords S.empty $
                     S.fromList ["type=experiment", "assay_term_name=" ++ x]

organismIs :: String -> KeyWords
organismIs x = KeyWords S.empty $
    S.fromList ["replicates.library.biosample.donor.organism.scientific_name=" ++ x]

cellIs :: String -> KeyWords
cellIs x = KeyWords S.empty $ S.fromList ["biosample_term_name=" ++ x]

-- * special search

openUrl :: String -> String -> IO B.ByteString
openUrl url datatype = do
    initReq <- parseUrl url
    let request = initReq { method = "GET" 
                          , requestHeaders = [("accept", "application/json")]
                          }
    r <- withManager $ \manager -> httpLbs request manager
    return $ responseBody r

jsonFromUrl :: String -> IO (Either String Value)
jsonFromUrl url = fmap eitherDecode $ openUrl (base ++ url) "application/json"

-- * Inspection

(|@) :: Value -> T.Text -> Value
(|@) (Object obj) key = M.lookupDefault (error "cannot find key") key obj
(|@) _ _ = error "not an object"
{-# INLINE (|@) #-}

(|!) :: Value -> Int -> Value
(|!) (Array ar) i = ar V.! i
(|!) _ _ = error "not an array"
{-# INLINE (|!) #-}

(&) :: a -> (a -> b) -> b
(&) = flip ($)
{-# INLINE (&) #-}

as :: FromJSON a => Value -> a
as = getResult . fromJSON
  where
    getResult (Error e) = error e 
    getResult (Success x) = x
{-# INLINE as #-}

test = do Right x <- search $ def <> assayIs "ChIP-seq" <> organismIs "Homo sapiens" <> cellIs "H1-hESC"
          showResult $ x !! 2
