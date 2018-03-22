{-# LANGUAGE OverloadedStrings #-}

module Bio.RealWorld.GENCODE
    ( Gene(..)
    , readGenes
    ) where

import           Conduit
import qualified Data.ByteString.Char8 as B
import           Data.CaseInsensitive  (CI, mk)
import qualified Data.HashMap.Strict   as M
import           Data.Maybe            (fromJust)

import           Bio.Utils.Misc        (readInt)

-- | GTF's position is 1-based, but here we convert it to 0-based indexing.
data Gene = Gene
    { geneName        :: !(CI B.ByteString)
    , geneId          :: !B.ByteString
    , geneChrom       :: !B.ByteString
    , geneLeft        :: !Int
    , geneRight       :: !Int
    , geneStrand      :: !Bool
    , geneTranscripts :: [(Int, Int)]
    } deriving (Show)

-- | Read gene information from Gencode GTF file
readGenes :: FilePath -> IO [Gene]
readGenes input = do
    (genes, transcripts) <- runResourceT $ runConduit $ sourceFile input .|
        linesUnboundedAsciiC .| foldlC f (M.empty, M.empty)
    return $ M.elems $ M.foldlWithKey'
        (\m k v -> M.adjust (\g -> g{geneTranscripts=v}) k m)
        genes transcripts
  where
    f (genes, transcripts) l
        | B.head l == '#' = (genes, transcripts)
        | f3 == "gene" =
            let g = Gene (mk $ getField "gene_name") i f1 (readInt f4 - 1)
                    (readInt f5 - 1) (f7=="+") []
                i = getField "gene_id"
            in (M.insert i g genes, transcripts)
        | f3 == "transcript" =
            let t = (readInt f4 - 1, readInt f5 - 1)
                i = getField "gene_id"
                updateFn Nothing = Just [t]
                updateFn (Just x) = Just $ t : x
            in (genes, M.alter updateFn i transcripts)
        | otherwise = (genes, transcripts)
      where
        [f1,_,f3,f4,f5,_,f7,_,f9] = B.split '\t' l
        fields = map (B.break (==' ') . strip) $ B.split ';' f9
        getField x = B.init $ B.drop 2 $ fromJust $ lookup x fields
    strip = fst . B.spanEnd isSpace . B.dropWhile isSpace
    isSpace = (== ' ')
