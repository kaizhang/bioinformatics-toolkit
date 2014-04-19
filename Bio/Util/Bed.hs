{-# LANGUAGE OverloadedStrings, UnicodeSyntax, BangPatterns #-}

module Bio.Util.Bed (
      readInt
    , readDouble
) where

import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString.Lex.Double as L
import Data.Maybe

readInt ∷ B.ByteString → Int
readInt = fst.fromJust.B.readInt

readDouble ∷ B.ByteString → Double
readDouble = fst . fromJust . L.readDouble
