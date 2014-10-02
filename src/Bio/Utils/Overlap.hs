module Bio.Utils.Overlap
    ( overlapFragment
    , overlapNucl
    , coverage
    ) where

import qualified Data.ByteString.Char8 as B
import qualified Data.IntervalMap.Strict as IM
import qualified Data.HashMap.Strict as M
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as VM
import Data.List
import Data.Function
import Bio.Data.Bed
import Control.Monad
import Data.Conduit
import qualified Data.Conduit.List as CL
import Control.Monad.Trans.Class (lift)

-- | convert lines of a BED file into a data structure - A hashmap of which the
-- | chromosomes, and values are interval maps.
toMap :: [(B.ByteString, (Int, Int))] -> M.HashMap B.ByteString (IM.IntervalMap Int Int)
toMap input = M.fromList.map create.groupBy ((==) `on` (fst.fst)) $ zip input [0..]
    where
        f ((_, x), i) = (toInterval x, i)
        create xs = (fst.fst.head $ xs, IM.fromDistinctAscList.map f $ xs)
{-# INLINE toMap #-}

{-
-- | coverages of bins
-- FIXME: Too ugly
coverage :: [BED]  -- ^ genomic locus in BED format
         -> [BED]  -- ^ reads in BED format
         -> (V.Vector Double, Int)
coverage bin tags = getResult (V.create (VM.replicate (n+1) 0 >>= go tags))
  where
    getResult v = (V.zipWith normalize (V.slice 0 n v) featWidth, v V.! n)
    go ts v = do
        forM_ ts (\t -> do
            let set = M.lookup (t^.chrom) featMap
                s = t^.chromStart
                e = t^.chromEnd
                b = (s, e)
                l = s - e + 1
                intervals = case set of
                    Just iMap -> IM.intersecting iMap.toInterval $ b
                    _ -> []
            forM_ intervals (\interval -> do 
                let i = snd interval
                    nucl = overlap b . fst $ interval
                VM.write v i . (+nucl) =<< VM.read v i
                )
            VM.write v n . (+l) =<< VM.read v n
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
    -}

coverage :: [BED]  -- ^ genomic locus in BED format
         -> Source IO BED  -- ^ reads in BED format
         -> IO (V.Vector Double, Int)
coverage bin tags = liftM getResult $ tags $$ sink
  where
    sink :: Sink BED IO (V.Vector Int)
    sink = do
        v <- lift $ VM.replicate (n+1) 0
        CL.mapM_ $ \t -> do
                let set = M.lookup (_chrom t) featMap
                    s = _chromStart t
                    e = _chromEnd t
                    b = (s, e)
                    l = e - s + 1
                    intervals = case set of
                        Just iMap -> IM.intersecting iMap.toInterval $ b
                        _ -> []
                forM_ intervals (\interval -> do 
                    let i = snd interval
                        nucl = overlap b . fst $ interval
                    VM.write v i . (+nucl) =<< VM.read v i
                    )
                VM.write v n . (+l) =<< VM.read v n
        lift $ V.freeze v
    getResult v = (V.zipWith normalize (V.slice 0 n v) featWidth, v V.! n)
    featMap = toMap.map (\x -> (_chrom x, (_chromStart x, _chromEnd x))) $ bin
    featWidth = V.fromList.map (\x -> _chromEnd x - _chromStart x) $ bin
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
toInterval (l, u) = IM.ClosedInterval l u
{-# INLINE toInterval #-}
