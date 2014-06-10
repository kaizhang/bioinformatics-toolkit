module Bio.Motif where

import Data.Array.Unboxed
import qualified Data.ByteString.Char8 as B
import Data.List

-- | position weight matrix for motifs
newtype PWM = PWM (UArray (Int, Int) Double)

score :: PWM -> B.ByteString -> Double
score (PWM pwm) sequ = foldl' (+) 0 . map f $ [0 .. len-1]
  where
    f i = case sequ `B.index` i of
              'A' -> logBase 10 (pwm ! (0, i))
              'C' -> logBase 10 (pwm ! (1, i))
              'G' -> logBase 10 (pwm ! (2, i))
              'T' -> logBase 10 (pwm ! (3, i))
              _   -> error "Bio.Motif.score: invalid nucleotide!"
    len = snd.snd.bounds $ pwm

