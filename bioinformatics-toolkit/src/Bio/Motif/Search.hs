{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}

module Bio.Motif.Search
    ( findTFBS
    , findTFBSWith
    , findTFBSSlow
    , maxMatchSc
    , optimalScoresSuffix
    ) where

import           Conduit
import           Control.Arrow               ((***))
import           Control.Monad.ST            (runST)
import           Control.Parallel.Strategies (parMap, rdeepseq)
import qualified Data.ByteString.Char8       as B
import           Data.Function               (on)
import qualified Data.HashMap.Strict         as HM
import qualified Data.HashSet                as S
import           Data.List                   (nubBy)
import qualified Data.Matrix.Unboxed         as M
import           Data.Ord                    (comparing)
import qualified Data.Vector.Unboxed         as U

import           Bio.Motif                   (Bkgd (..), Motif (..), PWM (..),
                                              pValueToScore, rcPWM, scores',
                                              size)
import           Bio.Seq                     (DNA, toBS)
import qualified Bio.Seq                     as Seq (length)

-- | given a user defined threshold, look for TF binding sites on a DNA
-- sequence, using look ahead search. This function doesn't search for binding
-- sites on the reverse strand
findTFBS :: Monad m
         => Bkgd
         -> PWM
         -> DNA a
         -> Double
         -> Bool    -- ^ whether to skip ambiguous sequences. Recommend: True
                    -- in most cases
         -> ConduitT i (Int, Double) m ()
findTFBS bg pwm dna thres skip = findTFBSWith sigma bg pwm dna thres skip
  where
    sigma = optimalScoresSuffix bg pwm
{-# INLINE findTFBS #-}

findTFBSWith :: Monad m
             => U.Vector Double   -- ^ best possible match score of suffixes
             -> Bkgd
             -> PWM
             -> DNA a
             -> Double
             -> Bool    -- ^ whether to skip ambiguous sequences. Recommend: True
                        -- in most cases
             -> ConduitT i (Int, Double) m ()
findTFBSWith sigma bg pwm dna thres skip = loop 0
  where
    loop !i | i >= l - n + 1 = return ()
            | otherwise = do let (d, sc) = searchFn bg pwm sigma dna i thres
                             if d == n - 1
                                then yield (i, sc) >> loop (i+1)
                                else loop (i+1)
    l = Seq.length dna
    n = size pwm
    searchFn | skip = lookAheadSearch'
             | otherwise = lookAheadSearch
{-# INLINE findTFBSWith #-}

{-
-- | a parallel version of findTFBS
findTFBS' :: Bkgd
          -> PWM
          -> DNA a
          -> Double
          -> Bool
          -> [Int]
findTFBS' bg pwm dna th skip = concat $ parMap rdeepseq f [0,step..l-n+1]
  where
    f x = loop x
      where
        loop i | i >= x+step || i >= l-n+1 = []
               | otherwise = let d = fst $ searchFn bg pwm sigma dna i th
                             in if d == n-1
                                   then i : loop (i+1)
                                   else loop (i+1)
    sigma = optimalScoresSuffix bg pwm
    l = Seq.length dna
    n = size pwm
    step = 500000
    searchFn | skip = lookAheadSearch'
             | otherwise = lookAheadSearch
{-# INLINE findTFBS' #-}
-}

-- | use naive search
findTFBSSlow :: Monad m => Bkgd -> PWM -> DNA a -> Double -> ConduitT i (Int, Double) m ()
findTFBSSlow bg pwm dna thres = scores' bg pwm dna .| loop 0
  where
    loop i = do v <- await
                case v of
                    Just v' ->  if v' >= thres then yield (i, v') >> loop (i+1)
                                               else loop (i+1)
                    _ -> return ()
{-# INLINE findTFBSSlow #-}

-- | the largest possible match scores starting from every position of a DNA sequence
maxMatchSc :: Bkgd -> PWM -> DNA a -> Double
maxMatchSc bg pwm dna = loop (-1/0) 0
  where
    loop !max' !i | i >= l - n + 1 = max'
                  | otherwise = if d == n - 1 then loop sc (i+1)
                                              else loop max' (i+1)
      where
        (d, sc) = lookAheadSearch bg pwm sigma dna i max'
    sigma = optimalScoresSuffix bg pwm
    l = Seq.length dna
    n = size pwm
{-# INLINE maxMatchSc #-}

optimalScoresSuffix :: Bkgd -> PWM -> U.Vector Double
optimalScoresSuffix (BG (a, c, g, t)) (PWM _ pwm) =
    U.fromList . tail . map (last sigma -) $ sigma
  where
    sigma = scanl f 0 $ M.toRows pwm
    f !acc xs = let (i, s) = U.maximumBy (comparing snd) .
                                U.zip (U.fromList ([0..3] :: [Int])) $ xs
                in acc + case i of
                    0 -> log $ s / a
                    1 -> log $ s / c
                    2 -> log $ s / g
                    3 -> log $ s / t
                    _ -> undefined
{-# INLINE optimalScoresSuffix #-}

lookAheadSearch :: Bkgd              -- ^ background nucleotide distribution
                -> PWM               -- ^ pwm
                -> U.Vector Double   -- ^ best possible match score of suffixes
                -> DNA a             -- ^ DNA sequence
                -> Int               -- ^ starting location on the DNA
                -> Double            -- ^ threshold
                -> (Int, Double)     -- ^ (d, sc_d), the largest d such that sc_d > th_d
lookAheadSearch (BG (a, c, g, t)) pwm sigma dna start thres = loop (0, -1) 0
  where
    loop (!acc, !th_d) !d
      | acc < th_d = (d-2, acc)
      | otherwise = if d >= n
                       then (d-1, acc)
                       else loop (acc + sc, thres - sigma U.! d) (d+1)
      where
        sc = case toBS dna `B.index` (start + d) of
              'A' -> log $! matchA / a
              'C' -> log $! matchC / c
              'G' -> log $! matchG / g
              'T' -> log $! matchT / t
              'N' -> 0
              'V' -> log $! (matchA + matchC + matchG) / (a + c + g)
              'H' -> log $! (matchA + matchC + matchT) / (a + c + t)
              'D' -> log $! (matchA + matchG + matchT) / (a + g + t)
              'B' -> log $! (matchC + matchG + matchT) / (c + g + t)
              'M' -> log $! (matchA + matchC) / (a + c)
              'K' -> log $! (matchG + matchT) / (g + t)
              'W' -> log $! (matchA + matchT) / (a + t)
              'S' -> log $! (matchC + matchG) / (c + g)
              'Y' -> log $! (matchC + matchT) / (c + t)
              'R' -> log $! (matchA + matchG) / (a + g)
              _   -> error "Bio.Motif.Search.lookAheadSearch: invalid nucleotide"
        matchA = addSome $ M.unsafeIndex (_mat pwm) (d,0)
        matchC = addSome $ M.unsafeIndex (_mat pwm) (d,1)
        matchG = addSome $ M.unsafeIndex (_mat pwm) (d,2)
        matchT = addSome $ M.unsafeIndex (_mat pwm) (d,3)
        addSome !x | x == 0 = pseudoCount
                   | otherwise = x
        pseudoCount = 0.0001
    n = size pwm
{-# INLINE lookAheadSearch #-}

-- | this version skip sequences contain ambiguous bases, like "N"
lookAheadSearch' :: Bkgd              -- ^ background nucleotide distribution
                 -> PWM               -- ^ pwm
                 -> U.Vector Double   -- ^ best possible match score of suffixes
                 -> DNA a             -- ^ DNA sequence
                 -> Int               -- ^ starting location on the DNA
                 -> Double            -- ^ threshold
                 -> (Int, Double)     -- ^ (d, sc_d), the largest d such that sc_d > th_d
lookAheadSearch' (BG (a, c, g, t)) pwm sigma dna start thres = loop (0, -1) 0
  where
    loop (!acc, !th_d) !d
      | acc < th_d = (d-2, acc)
      | otherwise = if d >= n
                       then (d-1, acc)
                       else loop (acc + sc, thres - sigma U.! d) (d+1)
      where
        sc = case toBS dna `B.index` (start + d) of
              'A' -> log $! matchA / a
              'C' -> log $! matchC / c
              'G' -> log $! matchG / g
              'T' -> log $! matchT / t
              _   -> -1 / 0
        matchA = addSome $ M.unsafeIndex (_mat pwm) (d,0)
        matchC = addSome $ M.unsafeIndex (_mat pwm) (d,1)
        matchG = addSome $ M.unsafeIndex (_mat pwm) (d,2)
        matchT = addSome $ M.unsafeIndex (_mat pwm) (d,3)
        addSome !x | x == 0 = pseudoCount
                   | otherwise = x
        pseudoCount = 0.0001
    n = size pwm
{-# INLINE lookAheadSearch' #-}

{-
data SpaceDistribution = SpaceDistribution
    { _motif1   :: Motif
    , _nSites1  :: (Int, Int)
    , _motif2   :: Motif
    , _nSites2  :: (Int, Int)
    , _same     :: [(Int, Int)]
    , _opposite :: [(Int, Int)]
    } deriving (Show, Read)

-- | search for spacing constraint between two TFs
spaceConstraint :: [(Motif, Motif)]   -- ^ motifs, names must be unique
                -> Bkgd    -- ^ backgroud nucleotide distribution
                -> Double  -- ^ p-Value for motif finding
                -> Int     -- ^ half window size, typical 5
                -> Int     -- ^ max distance to search, typical 300
                -> DNA a -> [SpaceDistribution]
spaceConstraint pairs bg th w k dna = flip map pairs $ \(a, b) ->
    let (m1, site1) = HM.lookupDefault undefined (_name a) sites
        (m2, site2) = HM.lookupDefault undefined (_name b) sites
        (same, opposite) = spaceConstraintHelper site1 site2 w k
    in SpaceDistribution m1 ((U.length *** U.length) site1) m2
        ((U.length *** U.length) site2) same opposite
  where
    motifs = nubBy ((==) `on` _name) $ concatMap (\(a,b) -> [a,b]) pairs
    findSites (Motif _ pwm) = (fwd, rev)
      where
        fwd = runST $ runConduit $ findTFBS bg pwm dna cutoff True .| sinkVector
        rev = runST $ runConduit $ findTFBS bg pwm' dna cutoff' True .|
              mapC (+ (size pwm - 1)) .| sinkVector
        cutoff = pValueToScore th bg pwm
        cutoff' = pValueToScore th bg pwm'
        pwm' = rcPWM pwm
    sites = HM.fromList $ map (\m -> (_name m, (m, findSites m))) motifs
{-# INLINE spaceConstraint #-}

spaceConstraintHelper :: (U.Vector Int, U.Vector Int)
                      -> (U.Vector Int, U.Vector Int)
                      -> Int
                      -> Int
                      -> ([(Int, Int)], [(Int, Int)])
spaceConstraintHelper (fw1, rv1) (fw2, rv2) w k = (same, oppose)
  where
    rs = let rs' = [-k, -k+2*w+1 .. 0]
         in rs' ++ map (*(-1)) (reverse rs')
    fw2' = S.fromList $ U.toList fw2
    rv2' = S.fromList $ U.toList rv2
    nOverlap :: U.Vector Int -> S.HashSet Int -> Int -> Int -> Int
    nOverlap xs ys w' i = U.foldl' f 0 xs
      where
        f acc x | any (`S.member` ys) [x + i - w' .. x + i + w'] = acc + 1
                | otherwise = acc
    same = zip rs $ zipWith (+) nFF nRR
      where
        nFF = map (nOverlap fw1 fw2' w) rs
        nRR = map (nOverlap rv1 rv2' w) $ reverse rs
    oppose = zip rs $ zipWith (+) nFR nRF
      where
        nFR = map (nOverlap fw1 rv2' w) rs
        nRF = map (nOverlap rv1 fw2' w) $ reverse rs
{-# INLINE spaceConstraintHelper #-}

computePValue :: Double -> [Int] -> [(Int, Double)]
computePValue p xs = zip xs $ map (pValue n p) xs
  where
    n = foldl' (+) 0 xs

pValue :: Int -> Double -> Int -> Double
pValue n p x | n > 2000 = complCumulative (poisson (fromIntegral n* p)) $ fromIntegral x
             | otherwise = complCumulative (binomial n p) $ fromIntegral x

-}