module Bio.HiC
    ( mkContactMap
    ) where

import Bio.SamTools.Bam
import qualified Bio.SamTools.BamIndex as BI
import Control.Monad (forM_)
import Control.Monad.IO.Class (liftIO)
import qualified Data.ByteString.Char8 as B
import Data.Conduit
import qualified Data.Conduit.List as CL
import Data.Maybe (fromJust)

import Bio.Data.Bam
import Bio.Data.Bed

-- | a specific base in the genome
type Site = (B.ByteString, Int)

-- | two sites interacting with each other form a contact
type Contact = (Site, Site)

-- | generate contanct map
mkContactMap :: FilePath -> BED3 -> Int -> Int -> Source IO ((Int, Int), Int)
mkContactMap bamFl (BED3 chr s e) w extend = do
    handle <- liftIO $ BI.open bamFl
    forM_ [0..n-1] $ \i ->
        forM_ [i..n-1] $ \j -> do
            let s1 = i * w
                c1 = s1 + w'
                s2 = j * w
                c2 = s2 + w'
            r <- liftIO $ readCount handle (w'+extend) ((chr, c1), (chr, c2))
            yield ((s1, s2), r)
  where
    n = (e - s) `div` w
    w' = w `div` 2

-- | the number of tag pairs that overlap with the region
-- covered by the contact
readCount :: BI.IdxHandle  -- ^ bam file handler
          -> Int           -- ^ half window width
          -> Contact
          -> IO Int
readCount handle w ((c1, p1), (c2, p2)) = viewBam handle (c1, p1-w, p1+w) $$ CL.fold f 0
  where
    f acc x =
        case mateTargetName x of
            Just chr -> if chr == c2
                           then let mp = fromIntegral . fromJust .  matePosition $ x
                                    l = fromIntegral . fromJust . mateTargetLen $ x
                                in if isOverlapped r2 (mp, mp+l)
                                      then acc + 1
                                      else acc
                           else acc
            _ -> acc
    isOverlapped (lo,hi) (lo',hi') = lo' < hi && hi' > lo
    r2 = (p2-w, p2+w)
{-# INLINE readCount #-}
