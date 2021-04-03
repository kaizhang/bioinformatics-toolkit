{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE LambdaCase #-}
module Bio.Data.Bam
    ( BAM
    , getBamHeader
    , streamBam
    , sinkBam
    , bamToBedC
    , bamToBed
    , bamToFragmentC
    , bamToFragment
    , bamToFastqC
    , bamToFastq
    , sortedBamToBedPE
    ) where

import Control.Monad (mzero)
import Data.List (foldl')
import           Bio.Data.Bed
import           Bio.Data.Fastq
import           Bio.HTS
import           Conduit

-- | Convert bam record to bed record. Unmapped reads will be discarded.
bamToBedC :: MonadIO m => BAMHeader -> ConduitT BAM BED m ()
bamToBedC header = mapC (bamToBed header) .| concatC
{-# INLINE bamToBedC #-}

-- | Convert bam record to fastq record.
bamToFastqC :: Monad m => ConduitT BAM Fastq m ()
bamToFastqC = mapC bamToFastq .| concatC
{-# INLINE bamToFastqC #-}

-- | Convert pairedend bam to fragment. 
bamToFragmentC :: Monad m => BAMHeader -> ConduitT BAM BED m ()
bamToFragmentC header = mapC (bamToFragment header) .| concatC
{-# INLINE bamToFragmentC #-}

bamToFragment :: BAMHeader -> BAM -> Maybe BED
bamToFragment header bam 
    | not (isFirstSegment flg) = Nothing
    | otherwise = do
        chr1 <- refName header bam
        chr2 <- mateRefName header bam
        if chr1 == chr2
            then return $ BED chr1 (min start1 start2) (max end1 end2)
                (Just $ queryName bam) Nothing Nothing
            else mzero
  where
    start1 = startLoc bam
    end1 = endLoc bam
    start2 = mateStartLoc bam
    end2 = mateStartLoc bam + ciglen cig
    cig = case queryAuxData ('M', 'C') bam of
        Just (AuxString x) -> string2Cigar x
        _ -> error "No MC tag. Please run samtools fixmate on file first."
    ciglen (CIGAR c) = foldl' f 0 c
        where f acc (n,x) = if x `elem` "MDN=X" then n + acc else acc
    flg = flag bam
{-# INLINE bamToFragment #-}

-- | Convert pairedend bam file to bed. the bam file must be sorted by names,
-- e.g., using "samtools sort -n". This condition is checked from Bam header.
sortedBamToBedPE :: Monad m => BAMHeader -> ConduitT BAM (BED, BED) m ()
sortedBamToBedPE header = case getSortOrder header of
    Queryname -> loopBedPE .| concatC
    _         -> error "Bam file must be sorted by NAME."
  where
    loopBedPE = (,) <$$> await <***> await >>= \case
        Nothing -> return ()
        Just (bam1, bam2) -> if queryName bam1 /= queryName bam2
            then error "Adjacent records have different query names. Aborted."
            else do
                yield $ (,) <$> bamToBed header bam1 <*> bamToBed header bam2
                loopBedPE
      where
        (<$$>) = fmap . fmap
        (<***>) = (<*>) . fmap (<*>)
{-# INLINE sortedBamToBedPE #-}

-- | Convert BAM to BED.
bamToBed :: BAMHeader -> BAM -> Maybe BED
bamToBed header bam = mkBed <$> refName header bam 
  where
    mkBed chr = BED chr start end nm sc str
    start = startLoc bam
    end = endLoc bam
    nm = Just $ queryName bam
    str = Just $ not $ isRev bam
    sc = Just $ fromIntegral $ mapq bam
{-# INLINE bamToBed #-}

-- | Convert BAM to Fastq.
bamToFastq :: BAM -> Maybe Fastq
bamToFastq bam = Fastq (queryName bam) <$> getSeq bam <*> qualityS bam
{-# INLINE bamToFastq #-}