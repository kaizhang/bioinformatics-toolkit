{-# LANGUAGE OverloadedStrings #-}
module Bio.GO.Parser
    ( readOWL
    ) where

import Control.Applicative ((<$>))
import qualified Data.ByteString.Lazy.Char8 as L
import Data.Maybe
import qualified Data.Text as T
import Data.Text.Encoding (encodeUtf8)
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
                  (encodeUtf8 . f $ findChild "oboInOwl:id" x)
                  (f $ findChild "oboInOwl:hasOBONamespace" x)
                  ((encodeUtf8 . T.replace "_" ":" . snd . T.breakOnEnd "/" . snd . head . getAttributes) <$> findChild "rdfs:subClassOf" x)
    f = getText . head . getChildren . fromJust
