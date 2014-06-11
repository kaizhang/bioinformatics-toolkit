module Bio.Motif 
    ( score
    , readPWM
    , scores
    ) where

import qualified Data.ByteString.Char8 as B
import Data.List
import Data.Default.Generics
import Bio.Utils.Bed (readDouble)
import Data.Matrix

-- | k x 4 position weight matrix for motifs
newtype PWM = PWM (Matrix Double)

newtype BkgdModel = BG (Double, Double, Double, Double)

instance Default BkgdModel where
    def = BG (0.25, 0.25, 0.25, 0.25)

-- | get scores of a long sequences at each position
scores :: BkgdModel -> PWM -> B.ByteString -> [Double]
scores bg p@(PWM pwm) = go
  where
    go s | B.length s >= len = score bg p (B.take len s) : go (B.tail s)
         | otherwise = []
    len = nrows pwm

score :: BkgdModel -> PWM -> B.ByteString -> Double
score (BG (a, c, g, t)) (PWM pwm) sequ = foldl' (+) 0 . map f $ [0 .. len-1]
  where
    f i = case sequ `B.index` i of
              'A' -> log $ (pwm ! (i, 0)) / a
              'C' -> log $ (pwm ! (i, 1)) / c
              'G' -> log $ (pwm ! (i, 2)) / g
              'T' -> log $ (pwm ! (i, 3)) / t
              _   -> error "Bio.Motif.score: invalid nucleotide!"
    len = nrows pwm

readPWM :: B.ByteString -> PWM
readPWM x = PWM $! fromLists . map (map readDouble.B.words) . B.lines $ x
