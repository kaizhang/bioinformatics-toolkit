{-# LANGUAGE OverloadedStrings, UnicodeSyntax #-}

import qualified Data.ByteString.Char8 as B
import System.Environment
import ChIP

run ∷ [FilePath] → IO ()
run [f1] = do
    tags ← B.readFile f1
    let x = [55..400]
    print $ naiveCC' tags x
run _ = error "incorrect arguments"

main ∷ IO ()
main = do
    args ← getArgs
    run args
