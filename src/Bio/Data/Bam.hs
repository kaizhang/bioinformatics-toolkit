--------------------------------------------------------------------------------
-- |
-- Module      :  $Header$
-- Copyright   :  (c) 2014 Kai Zhang
-- License     :  MIT

-- Maintainer  :  kai@kzhang.org
-- Stability   :  experimental
-- Portability :  portable

-- functions for processing BED files
--------------------------------------------------------------------------------

module Bio.Data.Bam
    ( readBam
    , bamToBed
    , viewBam
    ) where

import Bio.SamTools.Bam
import Bio.SamTools.BamIndex
import Control.Monad.Trans.Class (lift)
import qualified Data.ByteString.Char8 as B
import Data.Conduit (Source, Conduit, yield)
import qualified Data.Conduit.List as CL
import Data.Maybe (fromJust)
import Bio.Data.Bed

readBam :: FilePath -> Source IO Bam1
readBam fl = do handle <- lift $ openBamInFile fl
                go handle
  where
    go h = do x <- lift $ get1 h
              case x of
                  Nothing -> lift $ closeInHandle h
                  Just bam -> yield bam >> go h
{-# INLINE readBam #-}

bamToBed :: Monad m => Conduit Bam1 m BED
bamToBed = CL.mapMaybe f
  where
    f bam =case targetName bam of
        Just chr ->
            let start = fromIntegral . fromJust . position $ bam
                end = start + (fromIntegral . fromJust . queryLength) bam
                nm = Just . queryName $ bam
                strand = Just . not . isReverse $ bam
            in Just $ BED chr start end nm Nothing strand
        _ -> Nothing
{-# INLINE bamToBed #-}

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
