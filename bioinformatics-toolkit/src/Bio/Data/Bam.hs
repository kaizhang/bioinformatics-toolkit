{-# LANGUAGE FlexibleContexts #-}
module Bio.Data.Bam
    ( Bam
    , HeaderState
    , withBamFile
    , readBam
    , writeBam
    , bamToBedC
    , bamToBed
    , bamToFastqC
    , bamToFastq
    , sortedBamToBedPE
    ) where

import           Bio.Data.Bed
import           Bio.Data.Fastq
import           Bio.HTS
import           Bio.HTS.Types        (Bam, FileHeader (..))
import           Conduit
import           Control.Lens         ((&), (.~))
import           Control.Monad.Reader (ask, lift)

-- | Convert bam record to bed record. Unmapped reads will be discarded.
bamToBedC :: ConduitT Bam BED HeaderState ()
bamToBedC = mapMC bamToBed .| concatC
{-# INLINE bamToBedC #-}

-- | Convert bam record to fastq record.
bamToFastqC :: Monad m => ConduitT Bam Fastq m ()
bamToFastqC = mapC bamToFastq .| concatC
{-# INLINE bamToFastqC #-}

-- | Convert pairedend bam file to bed. the bam file must be sorted by names,
-- e.g., using "samtools sort -n". This condition is checked from Bam header.
sortedBamToBedPE :: ConduitT Bam (BED, BED) HeaderState ()
sortedBamToBedPE = do
    maybeBam <- await
    case maybeBam of
        Nothing -> return ()
        Just b' -> do
            leftover b'
            sortOrd <- getSortOrder <$> lift ask
            case sortOrd of
                Queryname -> loopBedPE .| concatC
                _         -> error "Bam file must be sorted by NAME."
  where
    loopBedPE = do
        pair <- (,) <$$> await <***> await
        case pair of
            Nothing -> return ()
            Just (bam1, bam2) -> if qName bam1 /= qName bam2
                then error "Adjacent records have different query names. Aborted."
                else do
                    bed1 <- lift $ bamToBed bam1
                    bed2 <- lift $ bamToBed bam2
                    yield $ (,) <$> bed1 <*> bed2
                    loopBedPE
      where
        (<$$>) = fmap . fmap
        (<***>) = (<*>) . fmap (<*>)
{-# INLINE sortedBamToBedPE #-}


bamToBed :: Bam -> HeaderState (Maybe BED)
bamToBed bam = do
    BamHeader hdr <- lift ask
    return $
        (\chr -> asBed chr start end & name .~ nm & score .~ sc & strand .~ str)
        <$> getChr hdr bam
  where
    start = fromIntegral $ position bam
    end = fromIntegral $ endPos bam
    nm = Just $ qName bam
    str = Just $ not $ isRev bam
    sc = Just $ fromIntegral $ mapq bam
{-# INLINE bamToBed #-}

bamToFastq :: Bam -> Maybe Fastq
bamToFastq bam = Fastq (qName bam) <$> getSeq bam <*> return mempty <*> qualityS bam
{-# INLINE bamToFastq #-}