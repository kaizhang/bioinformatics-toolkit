{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
module Bio.Motif 
    ( PWM (..)
    , subPWM
    , Motif (..)
    , readPWM
    , scores
    , score
    , findTFBS
    , toIUPAC
    , readMEME
    , writeFasta
    , readFasta

    -- * References
    -- $references
    ) where

import Prelude hiding (sum)
import qualified Data.ByteString.Char8 as B
import Data.Double.Conversion.ByteString
import qualified Data.Vector.Generic as G
import Data.Default.Generics
import Data.List (sortBy)
import Data.Ord (comparing)
import Bio.Utils.Bed (readDouble, readInt)
import Bio.Seq
import Numeric.LinearAlgebra.Data
import NLP.Scores
import Control.Monad.State.Lazy

-- | k x 4 position weight matrix for motifs
data PWM = PWM 
    { _nSites :: !(Maybe Int)  -- ^ number of sites used to generate this matrix
    , _mat :: !(Matrix Double)
    } deriving (Show)

-- | extract sub-PWM given starting position and length
subPWM :: Int -> Int -> PWM -> PWM
subPWM i l (PWM n mat) = PWM n $ subMatrix (i, 0) (l, 4) mat

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
        [a, b, _, _] = sortBy (flip (comparing snd)) $ zip "ACGT" $ G.toList v
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

-- | get scores of a long sequences at each position
scores :: BkgdModel -> PWM -> DNA a -> [Double]
scores bg p@(PWM _ pwm) dna = go $! toBS dna
  where
    go s | B.length s >= len = scoreHelp bg p (B.take len s) : go (B.tail s)
         | otherwise = []
    len = rows pwm

score :: BkgdModel -> PWM -> DNA a -> Double
score bg p dna = scoreHelp bg p $! toBS dna
{-# INLINE score #-}

-- | given a user defined threshold, look for TF binding sites on a DNA 
-- sequence. This function doesn't search for binding sites on the reverse strand
findTFBS :: Motif -> DNA a -> Double -> [Int]
findTFBS (Motif _ pwm) dna thres = go dna
  where
    go = fst . unzip . filter f . zip [1..] . scores def pwm
    f x = snd x >= thres * maxScore
    maxScore = sum . map (\x -> log (maximum x / 0.25)) . toLists . _mat $ pwm

scoreHelp :: BkgdModel -> PWM -> B.ByteString -> Double
scoreHelp (BG (a, c, g, t)) (PWM _ pwm) dna = sum . map f $ [0 .. len-1]
  where
    f i = case dna `B.index` i of
              'A' -> log $! matchA / a
              'C' -> log $! matchC / c
              'G' -> log $! matchG / g
              'T' -> log $! matchT / t
              'N' -> log $! (matchA / a + matchC / c + matchG / g + matchT / t) / 4
              'V' -> log $! (matchA / a + matchC / c + matchG / g) / 3
              'H' -> log $! (matchA / a + matchC / c + matchT / t) / 3
              'D' -> log $! (matchA / a + matchG / g + matchT / t) / 3
              'B' -> log $! (matchC / c + matchG / g + matchT / t) / 3
              'M' -> log $! (matchA / a + matchC / c) / 2
              'K' -> log $! (matchG / g + matchT / t) / 2
              'W' -> log $! (matchA / a + matchT / t) / 2
              'S' -> log $! (matchC / c + matchG / g) / 2
              'Y' -> log $! (matchC / c + matchT / t) / 2
              'R' -> log $! (matchA / a + matchG / g) / 2
              _   -> error "Bio.Motif.score: invalid nucleotide"
      where
        matchA = addSome $ pwm ! i ! 0
        matchC = addSome $ pwm ! i ! 1
        matchG = addSome $ pwm ! i ! 2
        matchT = addSome $ pwm ! i ! 3
        addSome !x | x == 0 = pseudoCount
                   | otherwise = x
        pseudoCount = 0.001
    len = rows pwm
{-# INLINE scoreHelp #-}

-- | read pwm from a matrix
readPWM :: B.ByteString -> PWM
readPWM x = PWM Nothing
              $ fromLists . map (map readDouble.B.words) . filter (not.B.null) . B.lines $ x

writePWM :: PWM -> B.ByteString
writePWM = B.unlines . map (B.unwords . map toShortest) . toLists . _mat

writeFasta :: FilePath -> [Motif] -> IO ()
writeFasta fl motifs = B.writeFile fl contents
  where
    contents = B.intercalate "" . map f $ motifs
    f x = B.unlines [">" `B.append` _name x, writePWM $ _pwm x]

readFasta :: FilePath -> IO [Motif]
readFasta fl = do contents <- B.readFile fl
                  return . map f . tail . B.split '>' $ contents
  where
    f x = let (nm, remain) = B.break (=='\n') x
          in Motif nm (readPWM remain)

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

-- $references
--
-- * Douglas R. Cavener. (1987) Comparison of the consensus sequence flanking
-- translational start sites in Drosophila and vertebrates.
-- /Nucleic Acids Research/ 15 (4): 1353–1361.
-- <doi:10.1093/nar/15.4.1353 http://nar.oxfordjournals.org/content/15/4/1353>
