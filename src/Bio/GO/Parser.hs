{-# LANGUAGE OverloadedStrings #-}
module Bio.GO.Parser
    ( readOWL
    , readOWLAsMap
    ) where

import Control.Arrow ((&&&))
import Control.Applicative ((<$>))
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.HashMap.Strict as M
import Data.Maybe
import qualified Data.Text as T
import Data.Text.Encoding (encodeUtf8)
import Text.XML.Expat.Proc
import Text.XML.Expat.Tree

import Bio.GO

readOWL :: FilePath -> IO [GO]
readOWL fl = do
    c <- L.readFile fl
    let (xml, _) = parse defaultParseOptions c
        goTerms = findChildren "owl:Class" (xml :: Node T.Text T.Text)
    return . map pickle $ goTerms
  where
    pickle x =
        let id' = encodeUtf8 . f $ findChild "oboInOwl:id" x
            label = f $ findChild "rdfs:label" x
            parent = ( encodeUtf8 . T.replace "_" ":" . snd
                     . T.breakOnEnd "/" . snd . head . getAttributes
                     ) <$> findChild "rdfs:subClassOf" x
            namespace = f $ findChild "oboInOwl:hasOBONamespace" x
        in GO id' label parent namespace
    f = getText . head . getChildren . fromJust

readOWLAsMap :: FilePath -> IO GOMap
readOWLAsMap fl = M.fromList . map (_oboId &&& id) <$> readOWL fl
