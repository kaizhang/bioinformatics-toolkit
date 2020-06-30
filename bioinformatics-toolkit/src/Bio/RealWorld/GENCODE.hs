{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}

module Bio.RealWorld.GENCODE
    ( Gene(..)
    , Transcript(..)
    , TranscriptType(..)
    , readGenes
    , readGenesC
    , getPromoters
    , getDomains
    ) where

import           Conduit
import qualified Data.ByteString.Char8 as B
import           Data.CaseInsensitive  (CI, mk)
import qualified Data.HashMap.Strict   as M
import Data.List.Ordered (nubSort)
import           Data.Maybe            (fromMaybe, fromJust, isNothing)
import Lens.Micro
import Data.List (foldl')
import Data.Char (toLower)
import qualified Data.Vector as V

import Bio.Data.Bed
import Bio.Data.Bed.Types
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
readGenes input = runResourceT $ runConduit $ sourceFile input .| readGenesC

readGenesC :: Monad m => ConduitT B.ByteString o m [Gene]
readGenesC = do
    (genes, transcripts, exons, utrs) <- readElements
    let t = M.fromList $ map (\(a,b) -> (transId b, (a,b))) transcripts
    return $ nubGene $ M.elems $ foldl' addTranscript
        (M.fromList $ map (\x -> (geneId x, x)) genes) $
        M.elems $ foldl' addUTR (foldl' addExon t exons) utrs
{-# INLINE readGenesC #-}

nubGene :: [Gene] -> [Gene]
nubGene gs = nubSort $ map nubG gs
  where
    nubG g = g { geneTranscripts = nubSort $ map nubT $ geneTranscripts g}
    nubT t = t { transExon = nubSort $ transExon t 
               , transUTR = nubSort $ transUTR t  }
{-# INLINE nubGene #-}

readElements :: Monad m => ConduitT B.ByteString o m
    ( [Gene], [(B.ByteString, Transcript)]
    , [(B.ByteString, (Int, Int))], [(B.ByteString, (Int, Int))] )
readElements = linesUnboundedAsciiC .| foldlC f ([], [], [], [])
  where
    f acc l
        | B.head l == '#' = acc
        | featType == "gene" = _1 %~ (gene:) $ acc
        | featType == "transcript" = _2 %~ ((gid, transcript):) $ acc
        | featType == "exon" = _3 %~ ((tid, exon):) $ acc
        | featType == "utr" = _4 %~ ((tid, utr):) $ acc
        | otherwise = acc
      where
        gene = Gene (mk $ fromMaybe (error "could not find \"gene_name\"") $
            getField "gene_name") gid chr lPos rPos (f7=="+") []
        transcript = Transcript tid lPos rPos (f7=="+") [] [] tTy
        exon = (lPos, rPos)
        utr = (lPos, rPos)
        [chr,_,f3,f4,f5,_,f7,_,f9] = B.split '\t' l
        gid = fromMaybe (error "could not find \"gene_id\"") $ getField "gene_id"
        tid = fromMaybe (error "could not find \"transcript_id\"") $ getField "transcript_id"
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

getPromoters :: Int   -- ^ upstream
             -> Int   -- ^ downstream
             -> Gene
             -> [BEDExt BED3 (Int, CI B.ByteString)]
getPromoters up down Gene{..} = map g $ nubSort tss
  where
    g x | geneStrand = BEDExt (asBed geneChrom (max 0 $ x - up) (x + down)) (x, geneName)
        | otherwise = BEDExt (asBed geneChrom (max 0 $ x - down) (x + up)) (x, geneName)
    tss | geneStrand = geneLeft : map transLeft geneTranscripts
        | otherwise = geneRight : map transRight geneTranscripts
{-# INLINE getPromoters #-}

-- | Compute genes' regulatory domains using the algorithm described in GREAT.
-- NOTE: the result doesn't contain promoters
getDomains :: BEDLike b
           => Int             -- ^ Extension length. A good default is 1M.
           -> [b] -- ^ A list of promoters
           -> [b] -- ^ Regulatory domains
getDomains ext genes
    | null genes = error "No gene available for domain assignment!"
    | otherwise = filter ((>0) . size) $ concatMap f $ triplet $
        [Nothing] ++ map Just basal ++ [Nothing]
  where
    f (left, Just bed, right) =
      [ chromStart .~ leftPos $ chromEnd .~ s $ bed
      , chromStart .~ e $ chromEnd .~ rightPos $ bed ]
      where
        chr = bed^.chrom
        s = bed^.chromStart
        e = bed^.chromEnd
        leftPos
            | isNothing left || chr /= fromJust left ^. chrom = max (s - ext) 0
            | otherwise = min s $ max (s - ext) $ fromJust left ^. chromEnd
        rightPos
            | isNothing right || chr /= fromJust right ^. chrom = e + ext   -- TODO: bound check
            | otherwise = max e $ min (e + ext) $ fromJust right ^. chromStart
    f _ = undefined
    triplet (x1:x2:x3:xs) = (x1,x2,x3) : triplet xs
    triplet _ = []
    basal = V.toList $ fromSorted $ sortBed genes
{-# INLINE getDomains #-}