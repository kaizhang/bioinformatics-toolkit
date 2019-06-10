{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE LambdaCase #-}
module Bio.Data.Bam
    ( BAM
    , getBamHeader
    , streamBam
    , sinkBam
    , bamToBedC
    , bamToBed
    , bamToFastqC
    , bamToFastq
    , sortedBamToBedPE
    ) where

import           Bio.Data.Bed
import           Bio.Data.Fastq
import           Bio.HTS
import           Bio.HTS.Types        (BAM, BAMHeader, SortOrder(..))
import           Conduit
import           Lens.Micro         ((&), (.~))

-- | Convert bam record to bed record. Unmapped reads will be discarded.
bamToBedC :: MonadIO m => BAMHeader -> ConduitT BAM BED m ()
bamToBedC header = mapC (bamToBed header) .| concatC
{-# INLINE bamToBedC #-}

-- | Convert bam record to fastq record.
bamToFastqC :: Monad m => ConduitT BAM Fastq m ()
bamToFastqC = mapC bamToFastq .| concatC
{-# INLINE bamToFastqC #-}

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
    mkBed chr = asBed chr start end &
        name .~ nm & score .~ sc & strand .~ str
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