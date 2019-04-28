{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}

module Bio.RealWorld.GENCODE
    ( Gene(..)
    , readGenes
    , streamElements
    ) where

import           Conduit
import qualified Data.ByteString.Char8 as B
import           Data.CaseInsensitive  (CI, mk)
import qualified Data.HashMap.Strict   as M
import Data.List.Ordered (nubSort)
import           Data.Maybe            (fromJust)

import Bio.Data.Bed.Types
import           Bio.Utils.Misc        (readInt)

-- | GTF's position is 1-based, but here we convert it to 0-based indexing.
data Gene = Gene
    { geneName        :: !(CI B.ByteString)
    , geneId          :: !B.ByteString
    , geneChrom       :: !B.ByteString
    , geneLeft        :: !Int
    , geneRight       :: !Int
    , geneStrand      :: !Bool
    , geneTranscripts :: ![(Int, Int)]
    , geneExon        :: ![(Int, Int)]
    } deriving (Show)

-- | Read gene information from Gencode GTF file
readGenes :: FilePath -> IO [Gene]
readGenes input = do
    (genes, transcripts, exon) <- runResourceT $ runConduit $ sourceFile input .|
        linesUnboundedAsciiC .| foldlC f (M.empty, M.empty, M.empty)
    return $ M.elems $ flip
        (M.foldlWithKey' (\m k v -> M.adjust (\g -> g{geneExon=nubSort v}) k m)) exon $
        flip (M.foldlWithKey' (\m k v -> M.adjust (\g -> g{geneTranscripts=nubSort v}) k m))
        transcripts genes
  where
    f (genes, transcripts, exon) l
        | B.head l == '#' = (genes, transcripts, exon)
        | f3 == "gene" =
            let g = Gene (mk $ getField "gene_name") gid f1 leftPos rightPos
                    (f7=="+") [] []
            in (M.insert gid g genes, transcripts, exon)
        | f3 == "transcript" = (genes, update transcripts, exon)
        | f3 == "exon" = (genes, transcripts, update exon)
        | otherwise = (genes, transcripts, exon)
      where
        gid = getField "gene_id"
        leftPos = readInt f4 - 1
        rightPos = readInt f5 - 1
        update m = let updateFn Nothing = Just [(leftPos, rightPos)]
                       updateFn (Just x) = Just $ (leftPos, rightPos) : x
                    in M.alter updateFn gid m
        [f1,_,f3,f4,f5,_,f7,_,f9] = B.split '\t' l
        fields = map (B.break isSpace . strip) $ B.split ';' f9
        getField x = B.init $ B.drop 2 $ fromJust $ lookup x fields
    strip = fst . B.spanEnd isSpace . B.dropWhile isSpace
    isSpace = (== ' ')

streamElements :: Monad m => ConduitT B.ByteString BED m ()
streamElements = linesUnboundedAsciiC .| concatMapC f
  where
    f l | B.head l == '#' = Nothing
        | otherwise = Just $ BED chr (readInt start - 1) (readInt end - 1)
            (Just name) Nothing (Just $ strand == "+")
      where
        [chr,_,name,start,end,_,strand,_,_] = B.split '\t' l