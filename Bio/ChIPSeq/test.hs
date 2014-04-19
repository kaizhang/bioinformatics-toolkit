{-# LANGUAGE OverloadedStrings, UnicodeSyntax #-}

import qualified Data.ByteString.Lazy.Char8 as B
import System.Environment
import qualified Data.HashSet as S
import ChIP
import Control.Concatenative
import qualified Data.Vector as V
import Control.Monad

run ∷ [FilePath] → IO ()
run [f1,f2] = do
    tags ← B.readFile f1
    maps ← B.readFile f2
    let (forwd, rev) = both S.fromList $ readBed $ B.lines tags
        m = V.fromList $ readInterval $ tail $ B.lines maps
    forM_ (S.toList forwd) (\i → do
        if binarySearch i m
           then return ()
           else print i
       )
run _ = error "incorrect arguments"

main ∷ IO ()
main = do
    args ← getArgs
    run args
