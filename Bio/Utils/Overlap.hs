module Bio.Utils.Overlap (
      overlapFragment
    , overlapNucl
    , coverage
) where

--import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.IntervalMap.Strict as IM
import qualified Data.HashMap.Strict as M
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as VM
import Data.List
import Data.Function
import Bio.Utils.Bed
import Control.Monad
import Control.Lens ((^.))

-- | convert lines of a BED file into a data structure - A hashmap of which the
-- | chromosomes, and values are interval maps.
toMap :: [(B.ByteString, (Int, Int))] -> M.HashMap B.ByteString (IM.IntervalMap Int Int)
{-# INLINE toMap #-}
toMap input = M.fromList.map create.groupBy ((==) `on` (fst.fst)) $ zip input [0..]
    where
        f ((_, x), i) = (toInterval x, i)
        create xs = (fst.fst.head $ xs, IM.fromDistinctAscList.map f $ xs)

-- | coverages of bins
coverage :: [BED]  -- ^ genomic locus in BED format
         -> [BED]  -- ^ reads in BED format
         -> V.Vector Double
coverage bin tags = V.zipWith normalize (V.create (VM.replicate n 0 >>= go tags)) featWidth
  where
    go ts v = do
        forM_ ts (\t -> do
            let set = M.lookup (t^.chrom) featMap
                b = (t^.chromStart, t^.chromEnd)
                intervals = case set of
                    Just iMap -> IM.intersecting iMap.toInterval $ b
                    _ -> []
            forM_ intervals (\interval -> do 
                let i = snd interval
                    nucl = overlap b . fst $ interval
                VM.write v i . (+nucl) =<< VM.read v i
                )
            )
        return v
    featMap = toMap.map (\x -> (x^.chrom, (x^.chromStart, x^.chromEnd))) $ bin
    featWidth = V.fromList.map (\x -> x^.chromEnd - x^.chromStart) $ bin
    n = length bin
    overlap (l, u) (IM.ClosedInterval l' u') 
        | l' >= l = if u' <= u then u'-l'+1 else u-l'+1
        | otherwise = if u' <= u then u'-l+1 else u-l+1
    overlap _ _ = 0
    normalize a b = fromIntegral a / fromIntegral b

overlapFragment, overlapNucl :: 
                  [(Int, Int)] -- ^ Ascending order list
                -> [(Int, Int)] -- ^ tags in any order
                -> V.Vector Int
overlapFragment xs ts = V.create (VM.replicate n 0 >>= go ts)
    where
        n = length xs
        iMap = IM.fromAscList $ zip (map toInterval xs) [0..]
        go ts' v = do
            forM_ ts' (\x -> do
                let indices = snd.unzip.IM.intersecting iMap.toInterval $ x
                forM_ indices (\i -> VM.write v i . (+1) =<< VM.read v i)
                )
            return v

overlapNucl xs ts = V.create (VM.replicate n 0 >>= go ts)
    where
        n = length xs
        iMap = IM.fromAscList $ zip (map toInterval xs) [0..]
        go ts' v = do
            forM_ ts' (\x -> do
                let intervals = IM.intersecting iMap.toInterval $ x
                forM_ intervals (\interval -> do 
                    let i = snd interval
                        nucl = overlap x . fst $ interval
                    VM.write v i . (+nucl) =<< VM.read v i
                    )
                )
            return v
        overlap (l, u) (IM.ClosedInterval l' u') 
            | l' >= l = if u' <= u then u'-l'+1 else u-l'+1
            | otherwise = if u' <= u then u'-l+1 else u-l+1
        overlap _ _ = 0

toInterval :: (a, a) -> IM.Interval a
{-# INLINE toInterval #-}
toInterval (l, u) = IM.ClosedInterval l u
