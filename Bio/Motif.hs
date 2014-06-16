{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE OverloadedStrings #-}
module Bio.Motif 
    ( score
    , readPWM
    , scores
    , toIUPAC
    , readMEME
    ) where

import Prelude hiding (sum)
import Data.Foldable (Foldable)
import qualified Data.ByteString.Char8 as B
import qualified Data.Vector as V
import Data.Default.Generics
import Bio.Utils.Bed (readDouble, readInt)
import Data.Matrix
import NLP.Scores
import Control.Monad.State.Lazy

-- | k x 4 position weight matrix for motifs
data PWM = PWM !(Maybe Int)  -- ^ number of sites
               !(Matrix Double)
    deriving (Show)

data Motif = Motif
    { _name :: !B.ByteString
    , _pwm :: !PWM
    } deriving (Show)

-- | background nucletide frequencies (A, C, G, T)
newtype BkgdModel = BG (Double, Double, Double, Double)

instance Default BkgdModel where
    def = BG (0.25, 0.25, 0.25, 0.25)

-- | convert pwm to consensus sequence
toIUPAC :: PWM -> String
toIUPAC (PWM _ pwm) = map f $ toRows pwm
  where
    f v | a == c && c == g && g == t    = 'N'
        | a == c && c == g && g == max' = 'V'
        | a == c && c == t && t == max' = 'H'
        | a == g && g == t && t == max' = 'D'
        | c == g && g == t && t == max' = 'B'
        | a == c && c == max'           = 'M'
        | g == t && t == max'           = 'K'
        | a == t && t == max'           = 'W'
        | c == g && g == max'           = 'S'
        | c == t && t == max'           = 'Y'
        | a == g && g == max'           = 'R'
        | a == max'                     = 'A'
        | c == max'                     = 'C'
        | g == max'                     = 'G'
        | t == max'                     = 'T'
        | otherwise                     = undefined
      where 
        [a, c, g, t] = V.toList v
        max' = V.maximum v

-- | calculate distance between PWMs
distanceBy :: (forall t. Foldable t => t Double -> t Double -> Double) -> PWM -> PWM -> Double
distanceBy f (PWM _ m1) (PWM _ m2) = sum $ zipWith f x y
  where (x, y) = (toRows m1, toRows m2)

-- Calculate distance by Kullback-Leibler divergence
deltaKL :: PWM -> PWM -> Double
deltaKL (PWM _ m1) (PWM _ m2) = sum $ zipWith kullbackLeibler x y
  where (x, y) = (toRows m1, toRows m2)

-- Calculate distance by Jensen-Shannon divergence
deltaJS :: PWM -> PWM -> Double
deltaJS (PWM _ m1) (PWM _ m2) = sum $ zipWith jensenShannon x y
  where (x, y) = (toRows m1, toRows m2)

-- | get scores of a long sequences at each position
scores :: BkgdModel -> PWM -> B.ByteString -> [Double]
scores bg p@(PWM _ pwm) = go
  where
    go s | B.length s >= len = score bg p (B.take len s) : go (B.tail s)
         | otherwise = []
    len = nrows pwm

score :: BkgdModel -> PWM -> B.ByteString -> Double
score (BG (a, c, g, t)) (PWM _ pwm) sequ = sum . map f $ [0 .. len-1]
  where
    f i = case sequ `B.index` i of
              'A' -> log $ addSome (pwm ! (i, 0)) / a
              'C' -> log $ addSome (pwm ! (i, 1)) / c
              'G' -> log $ addSome (pwm ! (i, 2)) / g
              'T' -> log $ addSome (pwm ! (i, 3)) / t
              _   -> error "Bio.Motif.score: invalid nucleotide!"
    len = nrows pwm
    addSome x | x == 0 = pseudoCount
              | otherwise = x
    pseudoCount = 0.001

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
               return $ if st == 2 then [toMotif str]
                                   else []
    startOfPwm = B.isPrefixOf "letter-probability matrix:"
    toMotif (name:n:xs) = Motif name pwm
      where
        pwm = PWM (Just $ readInt n) $ fromLists . map (map readDouble.B.words) $ xs
    toMotif _ = error "error"

toRows :: Matrix a -> [V.Vector a]
toRows m = map (`getRow` m) [1 .. nrows m]
