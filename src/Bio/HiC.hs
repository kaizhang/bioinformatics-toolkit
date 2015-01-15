module Bio.HiC
    ( mkContactMap
    , mkContactMap'
    ) where

import Bio.SamTools.Bam
import qualified Bio.SamTools.BamIndex as BI
import Control.Monad (forM_, when, liftM)
import Control.Monad.Primitive
import Control.Monad.Trans (lift)
import qualified Data.ByteString.Char8 as B
import Data.Conduit
import qualified Data.Conduit.List as CL
import Data.Maybe (fromJust)
import qualified Data.Vector as V
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Unboxed as U

import Bio.Data.Bam
import Bio.Data.Bed

-- | a specific base in the genome
type Site = (B.ByteString, Int)

-- | two sites interacting with each other form a contact
type Contact = (Site, Site)

{-
mkContactMap' :: FilePath -> [BED3] -> Int -> Int -> Source IO ..
mkContactMap' bamFl regions w extend = do
    handle <- liftIO $ BI.open bamFl
    forM_ [0..n-1] 
  where
    positions = V.scanl f 0 regions'
    n = V.last positions
    f acc (BED3 chr s e) = acc + (e - s) `div` w
    regions' = V.fromList regions
    g i = binarySearch positions i
    -}

-- | O(n * (logN + k)). n = number of bins, N = number of tags. Generate contanct
-- map using constant memory
mkContactMap :: FilePath -> BED3 -> Int -> Int -> Source IO ((Int, Int), Int)
mkContactMap bamFl (BED3 chr s e) w extend = do
    handle <- lift $ BI.open bamFl
    forM_ [0..n-1] $ \i ->
        forM_ [i..n-1] $ \j -> do
            let s1 = s + i * w
                c1 = s1 + w'
                s2 = s + j * w
                c2 = s2 + w'
            r <- lift $ readCount handle (w'+extend) ((chr, c1), (chr, c2))
            yield ((s1, s2), r)
  where
    n = (e - s) `div` w
    w' = w `div` 2

-- | O(N + n). Store matrix in memory, faster when region is big
mkContactMap' :: PrimMonad m => BED3 -> Int -> Int -> Sink Bam1 m (U.Vector Int)
mkContactMap' (BED3 chr s e) w extend = do
    vec <- lift $ GM.replicate vecLen 0
    loop vec
    lift . liftM (G.map (`div` 2)) . G.unsafeFreeze $ vec
  where
    n = (e - s) `div` w
    vecLen = n * (n+1) `div` 2
    loop v = do
        x <- await
        case x of
            Just bam ->
                let flag = do bamChr <- targetName bam
                              matChr <- mateTargetName bam
                              return $ bamChr == chr && matChr == chr
                in case flag of
                    Just True -> do
                        let (p, pMate) = getStarts bam
                            a = (p - s) `div` w
                            b = (pMate - s) `div` w
                            i | a < b = idx a b
                              | otherwise = idx b a
                        when (p >= s && pMate >= s && p < e && pMate < e) $ 
                            lift $ GM.read v i >>= GM.write v i . (+1)
                        loop v
                    _ -> loop v
            _ -> return ()
    idx i j | i == 0 = j
            | otherwise = i * (2 * n - i + 1) `div` 2 + j - i

-- | get starting location of bam and its mate
getStarts :: Bam1 -> (Int, Int)
getStarts bam = let p1 = fromIntegral . fromJust . position $ bam
                    p2 = fromIntegral . fromJust . matePosition $ bam
                    l = fromIntegral . fromJust . queryLength $ bam
                    p1' | isReverse bam = p1 + l
                        | otherwise = p1
                    p2' | isReverse bam = p2 + l
                        | otherwise = p2
                in (p1', p2')
{-# INLINE getStarts #-}

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
                                    l = fromIntegral . fromJust . queryLength $ x
                                in if isOverlapped r2 (mp, mp+l)
                                      then acc + 1
                                      else acc
                           else acc
            _ -> acc
    isOverlapped (lo,hi) (lo',hi') = lo' < hi && hi' > lo
    r2 = (p2-w, p2+w)
{-# INLINE readCount #-}
