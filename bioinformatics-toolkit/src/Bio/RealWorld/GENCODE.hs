{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}

module Bio.RealWorld.GENCODE
    ( Gene(..)
    , Transcript(..)
    , readGenes
    ) where

import           Conduit
import qualified Data.ByteString.Char8 as B
import           Data.CaseInsensitive  (CI, mk)
import qualified Data.HashMap.Strict   as M
import Data.List.Ordered (nubSort)
import           Data.Maybe            (fromJust)
import Lens.Micro
import Data.List (foldl')
import Data.Char (toLower)

import           Bio.Utils.Misc        (readInt)

data TranscriptType = Coding
                    | NonCoding
                    deriving (Show, Eq, Ord)

-- | GTF's position is 1-based, but here we convert it to 0-based indexing.
data Gene = Gene
    { geneName        :: !(CI B.ByteString)
    , geneId          :: !B.ByteString
    , geneChrom       :: !B.ByteString
    , geneLeft        :: !Int
    , geneRight       :: !Int
    , geneStrand      :: !Bool
    , geneTranscripts :: ![Transcript]
    } deriving (Show, Eq, Ord)

data Transcript = Transcript
    { transId          :: !B.ByteString
    , transLeft        :: !Int
    , transRight       :: !Int
    , transStrand      :: !Bool
    , transExon        :: ![(Int, Int)]
    , transUTR         :: ![(Int, Int)]
    , transType        :: TranscriptType
    } deriving (Show, Eq, Ord)

-- | Read gene information from Gencode GTF file
readGenes :: FilePath -> IO [Gene]
readGenes input = do
    (genes, transcripts, exons, utrs) <- readElements input
    let t = M.fromList $ map (\(a,b) -> (transId b, (a,b))) transcripts
    return $ nubGene $ M.elems $ foldl' addTranscript
        (M.fromList $ map (\x -> (geneId x, x)) genes) $
        M.elems $ foldl' addUTR (foldl' addExon t exons) utrs

nubGene :: [Gene] -> [Gene]
nubGene gs = nubSort $ map nubG gs
  where
    nubG g = g { geneTranscripts = nubSort $ map nubT $ geneTranscripts g}
    nubT t = t { transExon = nubSort $ transExon t 
               , transUTR = nubSort $ transUTR t  }
{-# INLINE nubGene #-}

readElements :: FilePath
             -> IO ( [Gene]
                   , [(B.ByteString, Transcript)]
                   , [(B.ByteString, (Int, Int))]
                   , [(B.ByteString, (Int, Int))] )
readElements input = runResourceT $ runConduit $ sourceFile input .|
    linesUnboundedAsciiC .| foldlC f ([], [], [], [])
  where
    f acc l
        | B.head l == '#' = acc
        | featType == "gene" = _1 %~ (gene:) $ acc
        | featType == "transcript" = _2 %~ ((gid, transcript):) $ acc
        | featType == "exon" = _3 %~ ((tid, exon):) $ acc
        | featType == "utr" = _4 %~ ((tid, utr):) $ acc
        | otherwise = acc
      where
        gene = Gene (mk $ fromJust $ getField "gene_name") gid chr lPos rPos (f7=="+") []
        transcript = Transcript tid lPos rPos (f7=="+") [] [] tTy
        exon = (lPos, rPos)
        utr = (lPos, rPos)
        [chr,_,f3,f4,f5,_,f7,_,f9] = B.split '\t' l
        gid = fromJust $ getField "gene_id"
        tid = fromJust $ getField "transcript_id"
        tTy = case getField "transcript_type" of
            Just "protein_coding" -> Coding
            Nothing -> Coding
            _ -> NonCoding
        lPos = readInt f4 - 1
        rPos = readInt f5 - 1
        featType = B.map toLower f3
        getField x = fmap (B.init . B.drop 2) $ lookup x $
            map (B.break isSpace . strip) $ B.split ';' f9
    strip = fst . B.spanEnd isSpace . B.dropWhile isSpace
    isSpace = (== ' ')
{-# INLINE readElements #-}

addExon :: M.HashMap B.ByteString (a, Transcript)
        -> (B.ByteString, (Int, Int))
        -> M.HashMap B.ByteString (a, Transcript)
addExon m (key, val) = M.adjust (\(x, trans) ->
    (x, trans{transExon = val : transExon trans})) key m
{-# INLINE addExon #-}

addUTR :: M.HashMap B.ByteString (a, Transcript)
       -> (B.ByteString, (Int, Int))
       -> M.HashMap B.ByteString (a, Transcript)
addUTR m (key, val) = M.adjust (\(x, trans) ->
    (x, trans{transUTR = val : transUTR trans})) key m
{-# INLINE addUTR #-}

addTranscript :: M.HashMap B.ByteString Gene
              -> (B.ByteString, Transcript)
              -> M.HashMap B.ByteString Gene
addTranscript m (key, val) = M.adjust (\gene ->
    gene{geneTranscripts = val : geneTranscripts gene}) key m
{-# INLINE addTranscript #-}

{-
streamElements :: Monad m => ConduitT B.ByteString BED m ()
streamElements = linesUnboundedAsciiC .| concatMapC f
  where
    f l | B.head l == '#' = Nothing
        | otherwise = Just $ BED chr (readInt start - 1) (readInt end - 1)
            (Just name) Nothing (Just $ strand == "+")
      where
        [chr,_,name,start,end,_,strand,_,_] = B.split '\t' l
-}