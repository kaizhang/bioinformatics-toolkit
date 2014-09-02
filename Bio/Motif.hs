{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
module Bio.Motif 
    ( PWM (..)
    , subPWM
    , rcPWM
    , Motif (..)
    , readPWM
    , scores
    , scores'
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
import Data.Conduit

-- | k x 4 position weight matrix for motifs
data PWM = PWM 
    { _nSites :: !(Maybe Int)  -- ^ number of sites used to generate this matrix
    , _mat :: !(Matrix Double)
    } deriving (Show)

-- | extract sub-PWM given starting position and length
subPWM :: Int -> Int -> PWM -> PWM
subPWM i l (PWM n mat) = PWM n $ subMatrix (i, 0) (l, 4) mat

-- | reverse complementary of PWM
rcPWM :: PWM -> PWM
rcPWM (PWM n m) = PWM n (fliprl . flipud $ m)

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

-- | a streaming version of scores
scores' :: Monad m => BkgdModel -> PWM -> DNA a -> Source m Double
scores' bg p@(PWM _ pwm) dna = go 0
  where
    go i | i < n - len + 1 = do yield $ scoreHelp bg p $ B.take len $ B.drop i s
                                go (i+1)
         | otherwise = return ()
    s = toBS dna
    n = B.length s
    len = rows pwm

score :: BkgdModel -> PWM -> DNA a -> Double
score bg p dna = scoreHelp bg p $! toBS dna
{-# INLINE score #-}

-- | given a user defined threshold (between 0 and 1), look for TF binding sites on a DNA 
-- sequence. This function doesn't search for binding sites on the reverse strand
findTFBS :: Monad m => Motif -> DNA a -> Double -> Source m Int
findTFBS (Motif _ pwm) dna thres = scores' def pwm dna
                                $= loop 0
  where
    loop i = do v <- await
                case v of
                    Just v' ->  if v' >= gate
                                   then yield i >> loop (i+1)
                                   else loop (i+1)
                    _ -> return ()
    gate = thres * maxScore
    maxScore = sum . map (\x -> log (maximum x / 0.25)) . toLists . _mat $ pwm
{-# INLINE findTFBS #-}

scoreHelp :: BkgdModel -> PWM -> B.ByteString -> Double
scoreHelp (BG (a, c, g, t)) (PWM _ pwm) dna = sum . map f $ [0 .. len-1]
  where
    f i = case dna `B.index` i of
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
      where
        matchA = addSome $ pwm ! i ! 0
        matchC = addSome $ pwm ! i ! 1
        matchG = addSome $ pwm ! i ! 2
        matchT = addSome $ pwm ! i ! 3
        addSome !x | x == 0 = pseudoCount
                   | otherwise = x
        pseudoCount = 0.0001
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
