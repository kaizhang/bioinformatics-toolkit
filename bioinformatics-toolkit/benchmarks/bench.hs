{-# LANGUAGE OverloadedStrings #-}

import Criterion.Main
import System.Random
import qualified Data.ByteString.Char8 as B
import Data.Default.Class
import Conduit
import Control.Monad.Identity
import System.IO.Unsafe
import AI.Clustering.Hierarchical
import Data.Either

import Bio.Data.Fasta
import Bio.Motif
import Bio.Motif.Search
import Bio.Motif.Alignment
import Bio.Seq
import Bio.Data.Fastq

dna :: DNA Basic
dna = fromRight undefined $ fromBS $ B.pack $ map f $ take 5000 $ randomRs (0, 3) (mkStdGen 2)
  where
    f :: Int -> Char
    f x = case x of
        0 -> 'A'
        1 -> 'C'
        2 -> 'G'
        3 -> 'T'
        _ -> undefined

pwm :: PWM
pwm = toPWM [ "0.3 0.3 0.3 0.1"
            , "0 0.5 0 0.5"
            , "0.1 0.2 0.5 0.3"
            , "0.1 0.1 0.1 0.7"
            , "0 0 0 1"
            , "0.25 0.25 0.25 0.25"
            , "0.1 0.1 0.3 0.5"
            , "0.25 0.25 0 0.5"
            , "0.1 0.1 0.7 0.1"
            , "0 0 0 1"
            ]

motifs :: [Motif]
motifs = unsafePerformIO $ readFasta' "data/motifs.fasta"

main :: IO ()
main = do
  fastqContent <- B.readFile "tests/data/test.fastq"
  defaultMain 
      [ bench "motif score" $ nf (scores def pwm) dna 
      --, bgroup "TFBS scanning" [ bench "Naive" $ nf (\x -> runIdentity $ findTFBSSlow def pwm x (0.6 * optimalScore def pwm) $$ CL.consume) dna
      --                         , bench "look ahead" $ nf (\x -> runIdentity $ findTFBS def pwm x (0.6 * optimalScore def pwm) $$ CL.consume) dna
      --                         ]
      --, bench "motif merge" $ nf (\x -> fmap show $ flatten $ buildTree x) motifs
      , bench "Read FASTQ" $ nfAppIO (\x -> runConduit $ yield x .| parseFastqC .| sinkList) fastqContent
      ]

