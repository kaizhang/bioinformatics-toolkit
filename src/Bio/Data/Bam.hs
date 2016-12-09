{-# LANGUAGE FlexibleContexts #-}
module Bio.Data.Bam
    ( Bam
    , HeaderState
    , runBam
    , readBam
    , writeBam
    , bamToBed
    , sortedBamToBedPE
    ) where

import           Bio.Data.Bed
import           Bio.HTS
import           Bio.HTS.Types             (Bam, FileHeader (..))
import           Conduit
import           Control.Monad.State (get, lift)

-- | Convert bam record to bed record. Unmapped reads will be discarded.
bamToBed :: Conduit Bam HeaderState BED
bamToBed = mapMC bamToBed1 =$= concatC
{-# INLINE bamToBed #-}

-- | Convert pairedend bam file to bed. the bam file must be sorted by names,
-- e.g., using "samtools sort -n". This condition is checked from Bam header.
sortedBamToBedPE :: Conduit Bam HeaderState (BED, BED)
sortedBamToBedPE = do
    maybeBam <- await
    case maybeBam of
        Nothing -> return ()
        Just b' -> do
            leftover b'
            sortOrd <- getSortOrder <$> lift get
            case sortOrd of
                Queryname -> loopBedPE =$= concatC
                _ -> error "Bam file must be sorted by NAME."
  where
    loopBedPE :: Conduit Bam HeaderState (Maybe (BED, BED))
    loopBedPE = do
        pair <- (,) <$$> await <***> await
        case pair of
            Nothing -> return ()
            Just (bam1, bam2) -> if qName bam1 /= qName bam2
                then error "Adjacent records have different query names. Aborted."
                else do
                    bed1 <- lift $ bamToBed1 bam1
                    bed2 <- lift $ bamToBed1 bam2
                    yield $ (,) <$> bed1 <*> bed2
                    loopBedPE
      where
        (<$$>) = fmap . fmap
        (<***>) = (<*>) . fmap (<*>)
{-# INLINE sortedBamToBedPE #-}


bamToBed1 :: Bam -> HeaderState (Maybe BED)
bamToBed1 bam = do
    BamHeader hdr <- lift get
    return $ (\chr -> BED chr start end nm sc strand) <$> getChr hdr bam
  where
    start = fromIntegral $ position bam
    end = fromIntegral $ endPos bam
    nm = Just $ qName bam
    strand = Just $ not $ isRev bam
    sc = Just $ fromIntegral $ mapq bam
{-# INLINE bamToBed1 #-}
