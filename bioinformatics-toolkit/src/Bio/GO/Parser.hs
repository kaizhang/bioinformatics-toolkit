{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE StrictData #-}
module Bio.GO.Parser
    ( readOWL
    , readOWLAsMap
    , GAF(..)
    , readGAF
    ) where

import           Control.Arrow              ((&&&))
import Conduit
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict        as M
import           Data.Text.Encoding         (decodeUtf8)
import           Text.XML.Expat.Proc
import           Text.XML.Expat.Tree
import qualified Data.CaseInsensitive as CI

import           Bio.GO
import           Bio.Utils.Misc (readInt)

readOWL :: FilePath -> IO [GO]
readOWL fl = do
    c <- L.readFile fl
    let (xml, _) = parse defaultParseOptions c
        goTerms = findChildren "owl:Class" (xml :: Node B.ByteString B.ByteString)
    return $ map process goTerms
  where
    process record = GO id' label parent namespace
      where
        id' = case findChild "oboInOwl:id" record of
            Nothing -> error "readOWL: cannot find id field"
            Just i -> readInt $ snd $ B.breakEnd (==':') $ getText $ head $ getChildren i
        label = case findChild "rdfs:label" record of
            Nothing -> error "readOWL: cannot find label field"
            Just l -> decodeUtf8 $ getText $ head $ getChildren l
        namespace = case findChild "oboInOwl:hasOBONamespace" record of
            Nothing -> error "readOWL: cannot find namespace field"
            Just ns -> decodeUtf8 $ getText $ head $ getChildren ns
        parent =
            let f p = case lookup "rdf:resource" (getAttributes p) of
                    Nothing -> error "readOWL: cannot find 'rdf:resource' attribute"
                    Just at -> readInt $ snd $ B.breakEnd (=='_') at
            in map f $ findChildren "rdfs:subClassOf" record 

readOWLAsMap :: FilePath -> IO GOMap
readOWLAsMap fl = M.fromListWith errMsg . map (_oboId &&& id) <$> readOWL fl
  where
    errMsg = error "readOWLAsMap: Duplicate records."


data GAF = GAF
    { gafDb :: B.ByteString
    , gafDbId :: B.ByteString
    , gafSymbol :: CI.CI B.ByteString
    , gafQualifier :: Maybe [B.ByteString]
    , gafGoId :: GOId
    , gafDbRef :: [B.ByteString]
    , gafEvidenceCode :: B.ByteString
    , gafWithOrFrom :: Maybe [B.ByteString]
    , gafAspect :: B.ByteString
    , gafName :: Maybe B.ByteString
    , gafSynonym :: Maybe [B.ByteString]
    , gafType :: B.ByteString
    , gafTaxon :: [B.ByteString]
    , gafDate :: B.ByteString
    , gafAssignedBy :: B.ByteString
    , gafAnnotationExtension :: Maybe [B.ByteString]
    , gafGeneProductID :: Maybe B.ByteString
    } deriving (Show, Ord, Eq)

-- | GO Annotation File (GAF) Format 2.1 Parser. For details read:
-- http://geneontology.org/page/go-annotation-file-gaf-format-21.
readGAF :: FilePath -> ConduitT i GAF (ResourceT IO) ()
readGAF input = sourceFileBS input .| linesUnboundedAsciiC .|
    (dropWhileC isCom >> mapC parseLine)
  where
    isCom l = B.head l == '!' || B.null l
{-# INLINE readGAF #-}

parseLine :: B.ByteString -> GAF
parseLine l = GAF f1 f2 (CI.mk f3) (optionals f4)
    (readInt $ snd $ B.breakEnd (==':') f5) (B.split '|' f6) f7 (optionals f8)
    f9 (optional f10) (optionals f11) f12 (B.split '|' f13) f14 f15
    (optionals f16) (optional f17)
  where
    [f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17] = B.split '\t' l
    optional x | B.null x = Nothing
               | otherwise = Just x
    optionals x | B.null x = Nothing
                | otherwise = Just $ B.split '|' x
{-# INLINE parseLine #-}
