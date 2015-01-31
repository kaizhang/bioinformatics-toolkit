{-# LANGUAGE BangPatterns #-}

module Bio.Motif.Search 
    ( findTFBS
    , findTFBS'
    , findTFBSSlow
    , maxMatchSc
    , MotifCompo(..)
    , spacingConstraint
    ) where

import Bio.Seq (DNA, toBS)
import qualified Bio.Seq as Seq (length)
import Bio.Motif
import Control.Parallel.Strategies (parMap, rdeepseq)
import Data.Conduit
import qualified Data.HashSet as S
import Data.List (foldl', intercalate)
import Data.Ord (comparing)
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector.Unboxed as U
import Statistics.Matrix hiding (map)


-- | given a user defined threshold, look for TF binding sites on a DNA 
-- sequence, using look ahead search. This function doesn't search for binding 
-- sites on the reverse strand
findTFBS :: Monad m => Bkgd -> PWM -> DNA a -> Double -> Source m Int
findTFBS bg pwm dna thres = loop 0
  where
    loop !i | i >= l - n + 1 = return ()
            | otherwise = do let (d, _) = lookAheadSearch bg pwm sigma dna i thres
                             if d == n - 1
                                then yield i >> loop (i+1)
                                else loop (i+1)
    sigma = optimalScoresSuffix bg pwm
    l = Seq.length dna
    n = size pwm
{-# INLINE findTFBS #-}

-- | a parallel version of findTFBS
findTFBS' :: Bkgd -> PWM -> DNA a -> Double -> [Int]
findTFBS' bg pwm dna th = concat $ parMap rdeepseq f [0,step..l-n+1]
  where
    f x = loop x
      where
        loop i | i >= x+step || i >= l-n+1 = []
               | otherwise = let d = fst $ lookAheadSearch bg pwm sigma dna i th
                             in if d == n-1
                                   then i : loop (i+1)
                                   else loop (i+1)
    sigma = optimalScoresSuffix bg pwm
    l = Seq.length dna
    n = size pwm
    step = 500000
{-# INLINE findTFBS' #-}

-- | use naive search
findTFBSSlow :: Monad m => Bkgd -> PWM -> DNA a -> Double -> Source m Int
findTFBSSlow bg pwm dna thres = scores' bg pwm dna $= 
                                       loop 0
  where
    loop i = do v <- await
                case v of
                    Just v' ->  if v' >= thres then yield i >> loop (i+1)
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
    sigma = scanl f 0 $ toRows pwm
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
        matchA = addSome $ unsafeIndex (_mat pwm) d 0
        matchC = addSome $ unsafeIndex (_mat pwm) d 1
        matchG = addSome $ unsafeIndex (_mat pwm) d 2
        matchT = addSome $ unsafeIndex (_mat pwm) d 3
        addSome !x | x == 0 = pseudoCount
                   | otherwise = x
        pseudoCount = 0.0001
    n = size pwm
{-# INLINE lookAheadSearch #-}

data MotifCompo = MotifCompo
    { _motif1 :: Motif
    , _motif2 :: Motif
    , _sameDirection :: Bool
    , _spacing :: Int
    , _score :: Double
    }

instance Show MotifCompo where
    show (MotifCompo m1 m2 isSame sp sc) = intercalate "\t"
        [B.unpack $ _name m1, B.unpack $ _name m2, show isSame, show sp, show sc]

-- | search for spacing constraint between two TFs
spacingConstraint :: Motif   -- ^ motif 1
                  -> Motif   -- ^ motif 2
                  -> Bkgd    -- ^ backgroud nucleotide distribution
                  -> Double  -- ^ p-Value for motif finding
                  -> Int     -- ^ half window size
                  -> Int     -- ^ max distance to search
                  -> DNA a -> [MotifCompo]
spacingConstraint m1@(Motif _ pwm1) m2@(Motif _ pwm2) bg th w k dna = same ++ oppose
  where
    rs = let rs' = [-k, -k+2*w+1 .. 0]
         in rs' ++ map (*(-1)) (reverse rs')
    -- on the same orientation
    same = zipWith f rs $ zipWith (+) nFF nRR
      where
        nFF = map (nOverlap forward1 forward2 w) rs
        nRR = map (nOverlap reverse1 reverse2 w) $ reverse rs
        f r c = MotifCompo m1 m2 True r $ fromIntegral c / n
    oppose = zipWith f rs $ zipWith (+) nFR nRF
      where
        nFR = map (nOverlap forward1 reverse2 w) rs
        nRF = map (nOverlap reverse1 forward2 w) $ reverse rs
        f r c = MotifCompo m1 m2 False r $ fromIntegral c / n
    forward1 = findTFBS' bg pwm1 dna th1
    forward2 = S.fromList $ findTFBS' bg pwm2 dna th2
    reverse1 = map (+ (s1 - 1)) $ findTFBS' bg (rcPWM pwm1) dna th1'
    reverse2 = S.fromList $ map (+ (s2 - 1)) $ findTFBS' bg (rcPWM pwm2) dna th2'
    th1 = pValueToScore' th bg pwm1
    th1' = pValueToScore' th bg $ rcPWM pwm1
    th2 = pValueToScore' th bg pwm2
    th2' = pValueToScore' th bg $ rcPWM pwm2
    s1 = size pwm1
    s2 = size pwm2
    n = sqrt $ fromIntegral $ (length forward1 + length reverse1) * (S.size forward2 + S.size reverse2)

    nOverlap :: [Int] -> S.HashSet Int -> Int -> Int -> Int
    nOverlap xs ys w' i = foldl' f 0 xs
      where
        f acc x | any (`S.member` ys) [x + i - w' .. x + i + w'] = acc + 1
                | otherwise = acc
{-# INLINE spacingConstraint #-}

{-
computePValue :: Double -> [Int] -> [(Int, Double)]
computePValue p xs = zip xs $ map (pValue n p) xs
  where
    n = foldl' (+) 0 xs
{-# INLINE computePValue #-}

pValue :: Int -> Double -> Int -> Double
pValue n p x | n > 2000 = complCumulative (poisson (fromIntegral n* p)) $ fromIntegral x
             | otherwise = complCumulative (binomial n p) $ fromIntegral x
{-# INLINE pValue #-}
-}
