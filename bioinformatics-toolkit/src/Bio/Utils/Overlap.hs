module Bio.Utils.Overlap
    ( overlapFragment
    , overlapNucl
    , coverage
    ) where

import           Bio.Data.Bed
import           Conduit
import           Lens.Micro                ((^.))
import           Control.Monad
import qualified Data.ByteString.Char8       as B
import           Data.Function
import qualified Data.HashMap.Strict         as M
import qualified Data.IntervalMap.Strict     as IM
import           Data.List
import qualified Data.Vector.Unboxed         as V
import qualified Data.Vector.Unboxed.Mutable as VM

-- | convert lines of a BED file into a data structure - A hashmap of which the
-- | chromosomes, and values are interval maps.
toMap :: [(B.ByteString, (Int, Int))] -> M.HashMap B.ByteString (IM.IntervalMap Int Int)
toMap input = M.fromList.map create.groupBy ((==) `on` (fst.fst)) $ zip input [0..]
    where
        f ((_, x), i) = (toInterval x, i)
        create xs = (fst.fst.head $ xs, IM.fromDistinctAscList.map f $ xs)
{-# INLINE toMap #-}

coverage :: [BED]  -- ^ genomic locus in BED format
         -> ConduitT () BED IO ()  -- ^ reads in BED format
         -> IO (V.Vector Double, Int)
coverage bin tags = liftM getResult $ runConduit $ tags .| sink
  where
    sink = do
        v <- lift $ VM.replicate (n+1) 0
        mapM_C $ \t -> do
                let set = M.lookup (t^.chrom) featMap
                    s = t^.chromStart
                    e = t^.chromEnd
                    b = (s, e)
                    l = e - s + 1
                    intervals = case set of
                        Just iMap -> IM.toList . IM.intersecting iMap . toInterval $ b
                        _ -> []
                forM_ intervals (\interval -> do
                    let i = snd interval
                        nucl = overlap b . fst $ interval
                    VM.write v i . (+nucl) =<< VM.read v i
                    )
                VM.write v n . (+l) =<< VM.read v n
        lift $ V.freeze v
    getResult v = (V.zipWith normalize (V.slice 0 n v) featWidth, v V.! n)
    featMap = toMap.map (\x -> (x^.chrom, (x^.chromStart, x^.chromEnd))) $ bin
    featWidth = V.fromList $ map size bin
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
                let indices = IM.elems . IM.intersecting iMap . toInterval $ x
                forM_ indices (\i -> VM.write v i . (+1) =<< VM.read v i)
                )
            return v

overlapNucl xs ts = V.create (VM.replicate n 0 >>= go ts)
    where
        n = length xs
        iMap = IM.fromAscList $ zip (map toInterval xs) [0..]
        go ts' v = do
            forM_ ts' (\x -> do
                let intervals = IM.toList . IM.intersecting iMap . toInterval $ x
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
toInterval (l, u) = IM.ClosedInterval l u
{-# INLINE toInterval #-}
