module Main where

import           Bio.Data.Bed
import           Bio.Data.Bed.Utils
import           Bio.Motif
import           Bio.Seq.IO
import           Conduit
import           Data.Default                      (def)
import           Data.Semigroup                    ((<>))
import           Data.Version                      (showVersion)
import           Options.Applicative
import           Paths_bioinformatics_toolkit_apps (version)
import Text.Printf
import System.IO (stdout)

data Options = Options
    { genomeFile :: FilePath
    , motifFile  :: FilePath
    , input      :: FilePath
    } deriving (Show, Read)

parser :: Parser Options
parser = Options
     <$> strArgument (metavar "GENOME")
     <*> strArgument (metavar "MOTIF_MEME")
     <*> strArgument
        ( metavar "INPUT"
       <> help "Motif binding sites in BED format, from running \"motifscan\".")

defaultMain :: Options -> IO ()
defaultMain opts = do
    withGenome (genomeFile opts) $ \genome -> do
        motifs <- readMEME $ motifFile opts
        runResourceT $ runConduit $ streamBed (input opts) .|
            getMotifScore genome motifs def .| sinkHandleBed stdout

main :: IO ()
main = execParser opts >>= defaultMain
  where
    opts = info (helper <*> parser) ( fullDesc <>
        header (printf "bioinformatics-toolkit-apps-v%s" (showVersion version)) )
