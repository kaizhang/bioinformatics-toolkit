{-# LANGUAGE BangPatterns      #-}
{-# LANGUAGE OverloadedStrings #-}

module Bio.Motif
    ( PWM(..)
    , size
    , subPWM
    , rcPWM
    , gcContentPWM
    , Motif(..)
    , Bkgd(..)
    , toPWM
    , ic
    , scores
    , scores'
    , score
    , optimalScore
    , CDF(..)
    , cdf
    , cdf'
    , truncateCDF
    , scoreCDF
    , pValueToScore
    , pValueToScoreExact
    , toIUPAC
    , readMEME
    , toMEME
    , fromMEME
    , writeMEME
    , writeFasta

    -- * References
    -- $references
    ) where

import           Conduit
import           Control.Arrow                     ((&&&))
import           Control.Monad.State.Lazy
import qualified Data.ByteString.Char8             as B
import           Data.Default.Class
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.List                         (foldl', sortBy)
import qualified Data.Matrix.Unboxed               as M
import           Data.Maybe                        (fromJust, isNothing)
import           Data.Ord                          (comparing)
import qualified Data.Vector.Algorithms.Intro      as I
import qualified Data.Vector.Unboxed               as U
import qualified Data.Vector.Unboxed.Mutable       as UM
import           Numeric.MathFunctions.Constants   (m_epsilon)
import           Prelude                           hiding (sum)
import           Text.Printf                       (printf)
import Statistics.Function (minMax)

import           Bio.Seq
import           Bio.Utils.Functions               (binarySearchBy)
import           Bio.Utils.Misc                    (readDouble, readInt)

-- | k x 4 position weight matrix for motifs
data PWM = PWM
    { _nSites :: !(Maybe Int)  -- ^ number of sites used to generate this matrix
    , _mat    :: !(M.Matrix Double)
    } deriving (Show, Read)

size :: PWM -> Int
size (PWM _ mat) = M.rows mat

-- | Information content of a poistion in pwm. (Not implemented)
ic :: PWM -> Int -> Double
ic = undefined

-- | Extract sub-PWM given starting position and length, zero indexed.
subPWM :: Int -> Int -> PWM -> PWM
subPWM i l (PWM n mat) = PWM n $ M.subMatrix (i,0) (i+l,3) mat
{-# INLINE subPWM #-}

-- | Reverse complementary of PWM.
rcPWM :: PWM -> PWM
rcPWM (PWM n mat) = PWM n . M.fromVector d . U.reverse . M.flatten $ mat
  where
    d = M.dim mat
{-# INLINE rcPWM #-}

-- | GC content of PWM.
gcContentPWM :: PWM -> Double
gcContentPWM (PWM _ mat) = loop 0 0 / fromIntegral m
  where
    m = M.rows mat
    loop !acc !i
        | i >= m = acc
        | otherwise =
            let acc' = acc + M.unsafeIndex mat (i,1) + M.unsafeIndex mat (i,2)
            in loop acc' (i+1)

data Motif = Motif
    { _name :: !B.ByteString
    , _pwm  :: !PWM
    } deriving (Show, Read)

-- | background model which consists of single nucletide frequencies, and di-nucletide
-- frequencies.
newtype Bkgd = BG (Double, Double, Double, Double)

instance Default Bkgd where
    def = BG (0.25, 0.25, 0.25, 0.25)

-- | Convert pwm to consensus sequence, see D. R. Cavener (1987).
toIUPAC :: PWM -> DNA IUPAC
toIUPAC (PWM _ pwm) = unsafeFromBS . B.pack . map f $ M.toRows pwm
  where
    f v | snd a > 0.5 && snd a > 2 * snd b = fst a
        | snd a + snd b > 0.75             = iupac (fst a, fst b)
        | otherwise                        = 'N'
      where
        [a, b, _, _] = sortBy (flip (comparing snd)) $ zip "ACGT" $ U.toList v
    iupac x = case sort' x of
        ('A', 'C') -> 'M'
        ('G', 'T') -> 'K'
        ('A', 'T') -> 'W'
        ('C', 'G') -> 'S'
        ('C', 'T') -> 'Y'
        ('A', 'G') -> 'R'
        _ -> undefined
    sort' (x, y) | x > y = (y, x)
                 | otherwise = (x, y)

-- | Get scores of a long sequences at each position.
scores :: Bkgd -> PWM -> DNA a -> [Double]
scores bg p@(PWM _ pwm) dna = go $! toBS dna
  where
    go s | B.length s >= len = scoreHelp bg p (B.take len s) : go (B.tail s)
         | otherwise = []
    len = M.rows pwm
{-# INLINE scores #-}

-- | A streaming version of scores.
scores' :: Monad m => Bkgd -> PWM -> DNA a -> ConduitT i Double m ()
scores' bg p@(PWM _ pwm) dna = go 0
  where
    go i | i < n - len + 1 = do yield $ scoreHelp bg p $ B.take len $ B.drop i s
                                go (i+1)
         | otherwise = return ()
    s = toBS dna
    n = B.length s
    len = M.rows pwm
{-# INLINE scores' #-}

score :: Bkgd -> PWM -> DNA a -> Double
score bg p dna = scoreHelp bg p $! toBS dna
{-# INLINE score #-}

-- | The best possible matching score of a pwm.
optimalScore :: Bkgd -> PWM -> Double
optimalScore (BG (a, c, g, t)) (PWM _ pwm) = foldl' (+) 0 . map f . M.toRows $ pwm
  where f xs = let (i, s) = U.maximumBy (comparing snd) .
                                U.zip (U.fromList ([0..3] :: [Int])) $ xs
               in case i of
                   0 -> log $ s / a
                   1 -> log $ s / c
                   2 -> log $ s / g
                   3 -> log $ s / t
                   _ -> undefined
{-# INLINE optimalScore #-}

scoreHelp :: Bkgd -> PWM -> B.ByteString -> Double
scoreHelp (BG (a, c, g, t)) (PWM _ pwm) dna = fst . B.foldl f (0,0) $ dna
  where
    f (!acc,!i) x = (acc + sc, i+1)
      where
        sc = case x of
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
              _   -> error "Bio.Motif.score: invalid nucleotide"
        matchA = addSome $ M.unsafeIndex pwm (i,0)
        matchC = addSome $ M.unsafeIndex pwm (i,1)
        matchG = addSome $ M.unsafeIndex pwm (i,2)
        matchT = addSome $ M.unsafeIndex pwm (i,3)
        addSome !y | y == 0 = pseudoCount
                   | otherwise = y
        pseudoCount = 0.0001
{-# INLINE scoreHelp #-}

-- | The probability that a kmer is generated by background based on a
-- 0-orderd Markov model.
pBkgd :: Bkgd -> B.ByteString -> Double
pBkgd (BG (a, c, g, t)) = B.foldl' f 1
  where
    f acc x = acc * sc
      where
        sc = case x of
              'A' -> a
              'C' -> c
              'G' -> g
              'T' -> t
              _ -> undefined
{-# INLINE pBkgd #-}

-- | calculate the minimum motif mathching score that would produce a kmer with
-- p-Value less than the given number. This score would then be used to search
-- for motif occurrences with significant p-Value
pValueToScore :: Double -> Bkgd -> PWM -> Double
pValueToScore p bg pwm = cdf' (scoreCDF bg pwm) $ 1 - p

-- | Unlike pValueToScore, this version compute the exact score but much slower
-- and is inpractical for long motifs.
pValueToScoreExact :: Double -- ^ desirable p-Value
              -> Bkgd
              -> PWM
              -> Double
pValueToScoreExact p bg pwm = go 0 0 $ sort' $
    map ((scoreHelp bg pwm &&& pBkgd bg) . B.pack) $ replicateM l "ACGT"
  where
    sort' xs = U.create $ do v <- U.unsafeThaw . U.fromList $ xs
                             I.sortBy (flip (comparing fst)) v
                             return v
    go !acc !i vec | acc > p = fst $ vec U.! (i - 1)
                   | otherwise = go (acc + snd (vec U.! i)) (i+1) vec
    l = size pwm

-- | The cumulative distribution function in the form of (x, P(X <= x)).
newtype CDF = CDF (U.Vector (Double, Double)) deriving (Read, Show)

-- P(X <= x)
cdf :: CDF -> Double -> Double
cdf (CDF v) x = let i = binarySearchBy cmp v x
                in case () of
                    _ | i >= n -> snd $ U.last v
                      | i == 0 -> if x == fst (U.head v) then snd (U.head v) else 0.0
                      | otherwise -> let (a, a') = v U.! (i-1)
                                         (b, b') = v U.! i
                                         α = (b - x) / (b - a)
                                         β = (x - a) / (b - a)
                                     in α * a' + β * b'
  where
    cmp (a,_) = compare a
    n = U.length v
{-# INLINE cdf #-}

-- | The inverse of cdf.
cdf' :: CDF -> Double -> Double
cdf' (CDF v) p
    | p > 1 || p < 0 = error "p must be in [0,1]"
    | otherwise =
        let i = binarySearchBy cmp v p
        in case () of
          _ | i >= n -> 1/0
            | i == 0 -> if p == snd (U.head v)
                then fst (U.head v)
                else error "cdf': impossible!"
            | otherwise -> let (a, a') = v U.! (i-1)
                               (b, b') = v U.! i
                               α = (b' - p) / (b' - a')
                               β = (p - a') / (b' - a')
                           in α * a + β * b
  where
    cmp (_,a) = compare a
    n = U.length v
{-# INLINE cdf' #-}

-- | Truncate the CDF by a value, in order to reduce the memory usage.
truncateCDF :: Double -> CDF -> CDF
truncateCDF x (CDF v) = CDF $ U.filter ((>=x) . snd) v
{-# INLINE truncateCDF #-}

-- | Approximate the cdf of motif matching scores using dynamic programming.
-- Algorithm:
-- Scan the PWM from left to right. For each position $i$, compute a score
-- density function $s_i$ such that $s_i(x)$ is the total number of sequences
-- with score $x$.
scoreCDF :: Bkgd -> PWM -> CDF
scoreCDF (BG (a,c,g,t)) pwm = toCDF $ loop (U.singleton 1, const 0) 0
  where
    loop (prev,scFn) i
        | i >= n = (prev, scFn)
        | lo < hi =
            let v = U.create $ do
                    new <- UM.replicate nBin' 0
                    flip U.imapM_ prev $ \x prob ->
                        when (prob /= 0) $ do
                            let sc = scFn x
                            UM.modify new (a * prob +) $ getIdx $ sc + a_at i
                            UM.modify new (c * prob +) $ getIdx $ sc + c_at i
                            UM.modify new (g * prob +) $ getIdx $ sc + g_at i
                            UM.modify new (t * prob +) $ getIdx $ sc + t_at i
                    return new
            in loop (v, \x -> (fromIntegral x + 0.5) * step + lo) (i+1)
        | otherwise = loop (prev, scFn) (i+1)
      where
        getIdx x = let j = truncate $ (x - lo) / step
                   in if j >= nBin' then nBin' - 1 else j
        lo = lo' + min'
        hi = hi' + max'
        nBin' = min 200000 $ ceiling $ (hi - lo) / precision
        step = (hi - lo) / fromIntegral nBin'
        lo' = scFn $ fst . U.head $ U.dropWhile ((==0) . snd) $ U.indexed prev
        hi' = scFn $ fst . U.head $ U.dropWhile ((==0) . snd) $ U.reverse $
            U.indexed prev
        (min', max') = minMax $ U.fromList [a_at i, c_at i, g_at i, t_at i]
    a_at i = log' (M.unsafeIndex (_mat pwm) (i,0)) - log a
    c_at i = log' (M.unsafeIndex (_mat pwm) (i,1)) - log c
    g_at i = log' (M.unsafeIndex (_mat pwm) (i,2)) - log g
    t_at i = log' (M.unsafeIndex (_mat pwm) (i,3)) - log t
    toCDF (v, scFn) = CDF $ compressCDF $ U.imap (\i x -> (scFn i, x)) $ U.scanl1 (+) v
    compressCDF v = U.ifilter f v
      where
        len = U.length v
        f i (_, x) | i == 0 || i == len = True
                   | otherwise = x - snd (v `U.unsafeIndex` (i-1)) > m_epsilon ||
                        snd (v `U.unsafeIndex` (i+1)) - x > m_epsilon
    precision = 1e-5
    n = size pwm
    log' x | x == 0 = log 0.001
           | otherwise = log x
{-# INLINE scoreCDF #-}

-- | Get pwm from a matrix.
toPWM :: [B.ByteString] -> PWM
toPWM x = PWM Nothing . M.fromLists . map (map readDouble.B.words) $ x

-- | Convert pwm to bytestring.
fromPWM :: PWM -> B.ByteString
fromPWM = B.unlines . map (B.unwords . map toShortest) . M.toLists . _mat

writeFasta :: FilePath -> [Motif] -> IO ()
writeFasta fl motifs = B.writeFile fl contents
  where
    contents = B.unlines . concatMap f $ motifs
    f x = [">" `B.append` _name x, fromPWM $ _pwm x]

readMEME :: FilePath -> IO [Motif]
readMEME = liftM fromMEME . B.readFile

writeMEME :: FilePath -> [Motif] -> Bkgd -> IO ()
writeMEME fl xs bg = B.writeFile fl $ toMEME xs bg

toMEME :: [Motif] -> Bkgd -> B.ByteString
toMEME xs (BG (a,c,g,t)) = B.intercalate "" $ header : map f xs
  where
    header = B.pack $ printf "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA %f C %f G %f T %f\n\n" a c g t
    f (Motif nm pwm) =
        let x = "MOTIF " `B.append` nm
            y = B.pack $ printf "letter-probability matrix: alength= 4 w= %d nsites= %d E= 0" (size pwm) sites
            z = B.unlines . map (B.unwords . ("":) . map toShortest) . M.toLists . _mat $ pwm
            sites | isNothing (_nSites pwm) = 1
                  | otherwise = fromJust $ _nSites pwm
        in B.unlines [x,y,z]
{-# INLINE toMEME #-}

fromMEME :: B.ByteString -> [Motif]
fromMEME meme = evalState (go $ B.lines meme) (0, [])
  where
    go :: [B.ByteString] -> State (Int, [B.ByteString]) [Motif]
    go (x:xs)
      | "MOTIF" `B.isPrefixOf` x = put (1, [B.drop 6 x]) >> go xs
      | otherwise = do
          (st, str) <- get
          case st of
              1 -> do when (startOfPwm x) $ put (2, str ++ [B.words x !! 7])
                      go xs
              2 -> let x' = B.dropWhile (== ' ') x
                   in if B.null x' || "URL" `B.isPrefixOf` x'
                         then do put (0, [])
                                 r <- go xs
                                 return (toMotif str : r)
                         else put (2, str ++ [x']) >> go xs
              _ -> go xs
    go [] = do (st, str) <- get
               return [toMotif str | st == 2]
    startOfPwm = B.isPrefixOf "letter-probability matrix:"
    toMotif (name:n:xs) = Motif name pwm
      where
        pwm = PWM (Just $ readInt n) $ M.fromLists . map (map readDouble.B.words) $ xs
    toMotif _ = error "error"
{-# INLINE fromMEME #-}

-- $references
--
-- * Douglas R. Cavener. (1987) Comparison of the consensus sequence flanking
-- translational start sites in Drosophila and vertebrates.
-- /Nucleic Acids Research/ 15 (4): 1353–1361.
-- <doi:10.1093/nar/15.4.1353 http://nar.oxfordjournals.org/content/15/4/1353>
