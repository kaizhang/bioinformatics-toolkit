{-# LANGUAGE FlexibleContexts #-}
module Bio.Data.Bam
    ( Bam
    , HeaderState
    , runBam
    , readBam
    , writeBam
    , bamToBed
    ) where

import           Bio.Data.Bed
import           Bio.HTS
import           Bio.HTS.Types             (Bam, FileHeader (..))
import           Conduit
import           Control.Monad.State (get, lift)

bamToBed :: Conduit Bam HeaderState BED
bamToBed = mapMC f =$= concatC
  where
    f bam = do
        BamHeader hdr <- lift get
        case getChr hdr bam of
            Just chr ->
                let start = fromIntegral $ position bam
                    end = fromIntegral $ endPos bam
                    nm = Just $ qName bam
                    strand = Just $ not $ isRev bam
                in return $ Just $ BED chr start end nm Nothing strand
            _ -> return Nothing
{-# INLINE bamToBed #-}

{-
viewBam :: IdxHandle -> (B.ByteString, Int, Int) -> Source IO Bam1
viewBam handle (chr, s, e) = case lookupTarget (idxHeader handle) chr of
    Nothing -> return ()
    Just chrId -> do
        q <- lift $ query handle chrId (fromIntegral s,fromIntegral e)
        go q
  where
    go q' = do r <- lift $ next q'
               case r of
                   Nothing -> return ()
                   Just bam -> yield bam >> go q'
{-# INLINE viewBam #-}
-}
