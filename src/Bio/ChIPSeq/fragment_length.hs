{-# LANGUAGE OverloadedStrings, UnicodeSyntax #-}

import System.Environment
import Bio.ChIPSeq.ChIP
import Bio.Utils.Bed

run ∷ [FilePath] → IO ()
run [f1] = do
    beds <- readBED f1
    let x = [55..400]
    print $ naiveCC' beds x
run _ = error "incorrect arguments"

main ∷ IO ()
main = do
    args ← getArgs
    run args
