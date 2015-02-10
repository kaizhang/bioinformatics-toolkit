{-# LANGUAGE OverloadedStrings #-}
module Bio.GO.Parser
    ( readOWL
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy.Char8 as L
import Data.Maybe
import qualified Data.Text as T
import Text.XML.Expat.Proc
import Text.XML.Expat.Tree

import Bio.GO

readOWL :: FilePath -> IO [GO]
readOWL fl = do c <- L.readFile fl
                let (xml, _) = parse defaultParseOptions c :: (Node T.Text T.Text, Maybe XMLParseError)
                    goTerms = findChildren "owl:Class" xml
                return $ map pickle goTerms
  where
    pickle x = GO (f $ findChild "rdfs:label" x)
                  (B.pack . T.unpack . f $ findChild "oboInOwl:id" x)
                  (f $ findChild "oboInOwl:hasOBONamespace" x)
    f = getText . head . getChildren . fromJust
