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
    ) where

import Bio.Data.Bed
import Bio.SamTools.Bam
import Control.Monad.Trans.Class
import Data.Conduit
import Data.Maybe

readBam :: FilePath -> Source IO Bam1
readBam fl = do handle <- lift $ openBamInFile fl
                go handle
  where
    go h = do x <- lift $ get1 h
              case x of
                  Nothing -> lift $ closeInHandle h
                  Just bam -> yield bam >> go h
{-# INLINE readBam #-}

bamToBed :: Conduit Bam1 IO BED
bamToBed = do
    x <- await
    case x of
        Nothing -> return ()
        Just bam -> case targetName bam of
            Nothing -> bamToBed
            Just chr -> do
                let start = fromIntegral . fromJust . position $ bam
                    end = start + (fromIntegral . fromJust . queryLength) bam - 1
                    nm = Just . queryName $ bam
                    strand = Just . not . isReverse $ bam
                yield $ BED chr start end nm Nothing strand
                bamToBed
{-# INLINE bamToBed #-}

