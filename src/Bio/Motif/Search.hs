{-# LANGUAGE BangPatterns #-}

module Bio.Motif.Search 
    ( findTFBS
    , findTFBS'
    ) where

import Bio.Seq
import Bio.Motif
import Data.Conduit
import Data.Ord
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector.Unboxed as U
import Statistics.Matrix hiding (map)

-- | given a user defined threshold, look for TF binding sites on a DNA 
-- sequence. This function doesn't search for binding sites on the reverse strand
findTFBS' :: Monad m => Bkgd -> Motif -> DNA a -> Double -> Source m Int
findTFBS' bg (Motif _ pwm) dna thres = scores' bg pwm dna $= 
                                       loop 0
  where
    loop i = do v <- await
                case v of
                    Just v' ->  if v' >= thres
                                   then yield i >> loop (i+1)
                                   else loop (i+1)
                    _ -> return ()
{-# INLINE findTFBS' #-}

findTFBS :: Monad m => Bkgd -> Motif -> DNA a -> Double -> Source m Int
findTFBS bg@(BG (a, c, g, t)) (Motif _ pwm) dna thres = loop 0
  where
    loop !i | i >= l - n + 1 = return ()
            | otherwise = do let (d, _) = lookAheadSearch bg pwm sigma dna i thres
                             if d == n - 1
                                then yield i >> loop (i+1)
                                else loop (i+1)
    sigma = U.fromList . tail . map ((-) (last sigma')) $ sigma'
    sigma' = scanl f 0 $ toRows $ _mat pwm
    l = Bio.Seq.length dna
    n = size pwm
    f !acc xs = let (i, s) = U.maximumBy (comparing snd) . 
                                U.zip (U.fromList ([0..3] :: [Int])) $ xs
                in acc + case i of
                    0 -> log $ s / a
                    1 -> log $ s / c
                    2 -> log $ s / g
                    3 -> log $ s / t
                    _ -> undefined
{-# INLINE findTFBS #-}

lookAheadSearch :: Bkgd -> PWM -> U.Vector Double -> DNA a -> Int -> Double -> (Int, Double)
lookAheadSearch (BG (a, c, g, t)) pwm sigma dna start thres = loop (0, -1) 0
  where
    loop (!acc, !th_d) !d 
      | acc < th_d = (d-2, acc)
      | otherwise = if d >= n
                       then (d-1, acc)
                       else loop (acc + sc, thres - sigma U.! d) (d+1)
      where
        sc = case (toBS dna) `B.index` (start + d) of
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

toRows :: Matrix -> [Vector]
toRows (Matrix _ ncol _ v) = loop v 
  where 
    loop x | U.length x >= ncol = a : loop b
           | otherwise = []
      where (a, b) = U.splitAt ncol x
{-# INLINE toRows #-}
