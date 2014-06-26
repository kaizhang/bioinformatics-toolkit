{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
module Bio.Motif 
    ( PWM (..)
    , subPWM
    , deltaKL
    , deltaJS
    , Motif (..)
    , readPWM
    , scores
    , score
    , toIUPAC
    , readMEME

    -- * References
    -- $references
    ) where

import Prelude hiding (sum)
import Data.Traversable (Traversable)
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector as V
import Data.Default.Generics
import Data.List (sortBy)
import Data.Ord (comparing)
import Bio.Utils.Bed (readDouble, readInt)
import Bio.Seq
import Data.Matrix
import NLP.Scores
import Control.Monad.State.Lazy

-- | k x 4 position weight matrix for motifs
data PWM = PWM 
    { _nSites :: !(Maybe Int)  -- ^ number of sites used to generate this matrix
    , _mat :: !(Matrix Double)
    } deriving (Show)

-- | extract sub-PWM given starting position and length, 1-based indexing
subPWM :: Int -> Int -> PWM -> PWM
subPWM i l (PWM n mat) = PWM n $ submatrix i (i+l-1) 1 4 mat

data Motif = Motif
    { _name :: !B.ByteString
    , _pwm :: !PWM
    } deriving (Show)

-- | background nucletide frequencies (A, C, G, T)
newtype BkgdModel = BG (Double, Double, Double, Double)

instance Default BkgdModel where
    def = BG (0.25, 0.25, 0.25, 0.25)

-- | convert pwm to consensus sequence, see D. R. Cavener (1987).
toIUPAC :: PWM -> DNA IUPAC
toIUPAC (PWM _ pwm) = fromBS . B.pack . map f $ toRows pwm
  where
    f v | snd a > 0.5 && snd a > 2 * snd b = fst a
        | snd a + snd b > 0.75             = iupac (fst a, fst b)
        | otherwise                        = 'N'
      where 
        [a, b, _, _] = sortBy (flip (comparing snd)) $ zip "ACGT" $ V.toList v
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

-- | calculate distance between PWMs
distanceBy :: (Traversable t, t ~ V.Vector) => (t Double -> t Double -> Double) -> PWM -> PWM -> Double
distanceBy f (PWM _ m1) (PWM _ m2) = mean $ zipWith f x y
  where (x, y) = (toRows m1, toRows m2)
{-# INLINE distanceBy #-}

-- Calculate distance by Kullback-Leibler divergence
deltaKL :: PWM -> PWM -> Double
deltaKL = distanceBy kullbackLeibler
{-# INLINE deltaKL #-}

-- Calculate distance by Jensen-Shannon divergence
deltaJS :: PWM -> PWM -> Double
deltaJS = distanceBy jensenShannon
{-# INLINE deltaJS #-}

-- | get scores of a long sequences at each position
scores :: BkgdModel -> PWM -> DNA a -> [Double]
scores bg p@(PWM _ pwm) dna = go $! toBS dna
  where
    go s | B.length s >= len = scoreHelp bg p (B.take len s) : go (B.tail s)
         | otherwise = []
    len = nrows pwm

score :: BkgdModel -> PWM -> DNA a -> Double
score bg p dna = scoreHelp bg p $ toBS dna

scoreHelp :: BkgdModel -> PWM -> B.ByteString -> Double
scoreHelp (BG (a, c, g, t)) (PWM _ pwm) dna = sum . map f $ [0 .. len-1]
  where
    f i = case dna `B.index` i of
              'A' -> log $ matchA / a
              'C' -> log $ matchC / c
              'G' -> log $ matchG / g
              'T' -> log $ matchT / t
              'N' -> log $ (matchA / a + matchC / c + matchG / g + matchT / t) / 4
              'V' -> log $ (matchA / a + matchC / c + matchG / g) / 3
              'H' -> log $ (matchA / a + matchC / c + matchT / t) / 3
              'D' -> log $ (matchA / a + matchG / g + matchT / t) / 3
              'B' -> log $ (matchC / c + matchG / g + matchT / t) / 3
              'M' -> log $ (matchA / a + matchC / c) / 2
              'K' -> log $ (matchG / g + matchT / t) / 2
              'W' -> log $ (matchA / a + matchT / t) / 2
              'S' -> log $ (matchC / c + matchG / g) / 2
              'Y' -> log $ (matchC / c + matchT / t) / 2
              'R' -> log $ (matchA / a + matchG / g) / 2
              _   -> error "Bio.Motif.score: invalid nucleotide"
      where
        matchA = addSome $! pwm ! (i+1, 1)
        matchC = addSome $! pwm ! (i+1, 2)
        matchG = addSome $! pwm ! (i+1, 3)
        matchT = addSome $! pwm ! (i+1, 4)
        addSome x | x == 0 = pseudoCount
                  | otherwise = x
        pseudoCount = 0.001
    len = nrows pwm
{-# INLINE scoreHelp #-}

-- | read pwm from a matrix
readPWM :: B.ByteString -> PWM
readPWM x = PWM Nothing
              $ fromLists . map (map readDouble.B.words) . B.lines $ x

readMEME :: FilePath -> IO [Motif]
readMEME = liftM fromMEME . B.readFile

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
                   in if B.null x'
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
        pwm = PWM (Just $ readInt n) $ fromLists . map (map readDouble.B.words) $ xs
    toMotif _ = error "error"
{-# INLINE fromMEME #-}

toRows :: Matrix a -> [V.Vector a]
toRows m = map (`getRow` m) [1 .. nrows m]

-- $references
--
-- * Douglas R. Cavener. (1987) Comparison of the consensus sequence flanking
-- translational start sites in Drosophila and vertebrates.
-- /Nucleic Acids Research/ 15 (4): 1353–1361.
-- <doi:10.1093/nar/15.4.1353 http://nar.oxfordjournals.org/content/15/4/1353>
