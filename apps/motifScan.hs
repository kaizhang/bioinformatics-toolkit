{-# LANGUAGE OverloadedStrings #-}
import Bio.Data.Fasta
import Bio.Seq
import Bio.Motif.Search
import Bio.Motif
import Data.Maybe
import Conduit
import Control.Monad.Identity
import System.Environment
import qualified Data.ByteString.Char8 as B

main = do
    [motifFl, dnaFl, motifId, p, a, c, g, t] <- getArgs
    motifs <- readMEME motifFl
    let pwm = fromJust $ lookup (B.pack motifId) $ map (\x -> (_name x, _pwm x)) motifs
        pwm' = rcPWM pwm
        cutoff = pValueToScore (read p) bg pwm
        cutoff' = pValueToScore (read p) bg' pwm'
        bg = BG (read a, read c, read g, read t)
        bg' = BG (read t, read g, read c, read a)
        f (name, xs) = let dna = fromBS (mconcat xs) :: DNA IUPAC
                       in name `B.append` "\t" `B.append`
                          (B.intercalate "," $ map (B.pack . show) $ runIdentity $
                          (findTFBS bg pwm dna cutoff True =$= mapC (+1)) >>
                          (findTFBS bg' pwm' dna cutoff' True =$= mapC (negate . (+1))) $$
                          sinkList)
    fastaReader dnaFl $$ mapM_C (B.putStrLn . f)
