{-# LANGUAGE OverloadedStrings, UnicodeSyntax, BangPatterns #-}

import Bio.Util.Overlap
import Bio.Util.Bed
import Criterion.Main

intervals ∷ [(Int,Int)]
intervals = binBySize (1, 100000) 500

tags ∷ [(Int,Int)]
tags = binBySize (1,1000000) 10

main ∷ IO ()
main = defaultMain [
      bench "bench1" $ whnf (overlapFragment intervals) tags
    ]
