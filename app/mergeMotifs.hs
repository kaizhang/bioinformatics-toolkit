{-# LANGUAGE OverloadedStrings #-}

import           AI.Clustering.Hierarchical hiding (drawDendrogram)
import           Control.Monad              (forM)
import qualified Data.ByteString.Char8      as B
import           Data.Default.Class
import           Data.List
import           Data.List.Split            (splitOn)
import           Data.Ord
{-
import           Diagrams.Backend.Cairo
import           Diagrams.Plots.Dendrogram
import           Diagrams.Prelude           (dims2D, strutX, (|||))
-}
import           Options.Applicative
import           System.IO
import           Text.Printf

import           Bio.Data.Fasta
import           Bio.Motif
import           Bio.Motif.Alignment
import           Bio.Motif.Merge
import           Bio.Seq                    (toBS)
import           Bio.Utils.Functions

data Options = Options
    { input    :: FilePath
    , output   :: FilePath
    , prefix   :: String
    , thres    :: Double
    , mode     :: String
    , svg      :: Maybe FilePath
    , dumpDist :: Bool
    , gap      :: Double
    } deriving (Show)

parser :: Parser Options
parser = Options
     <$> argument str (metavar "INPUT")
     <*> strOption
           ( long "output"
          <> short 'o'
          <> value "merged_output.meme"
          <> metavar "OUTPUT" )
     <*> strOption
           ( long "prefix"
          <> short 'p'
          <> value "merged"
          <> metavar "PREFIX"
          <> help "PREFIX that being add to the name of the merged motif" )
     <*> option auto
           ( long "thres"
          <> short 't'
          <> value 0.2
          <> metavar "THRESHOLD"
          <> help "two motifs that have distance belowing threshold would be merged, default is 0.2" )
     <*> strOption
           ( long "mode"
          <> short 'm'
          <> value "iter"
          <> metavar "MODE"
          <> help "merging algorithm, could be iter or tree, default is iter" )
     <*> (optional . strOption)
           ( long "svg"
          <> metavar "SVG"
          <> help "draw merging tree in svg format, only available in tree mode" )
     <*> switch
           ( long "dump-dist"
          <> help "output pairwise distances of original motifs without performing any merging" )
     <*> option auto
           ( long "gap"
          <> short 'g'
          <> value 0.15
          <> metavar "GAP_PENALTY"
          <> help "gap penalty" )

-- | Return pairwise distances
pairDistance :: [Motif] -> [(B.ByteString, B.ByteString, Double)]
pairDistance ms = map (\(a,b) -> (_name a, _name b, fst $ alignment (_pwm a) (_pwm b))) $ comb ms

treeMerge :: Double -> String -> [Motif] -> Double -> ([Motif], Dendrogram Motif)
treeMerge th pre ms gap = (zipWith f [0::Int ..] $ map merge $ tree `cutAt` th, tree)
  where
    f i (suffix, pwm) = Motif ((B.pack $ pre ++ "_" ++ show i ++ "_" ++ show (toIUPAC pwm))
                                 `B.append` "("
                                 `B.append` suffix
                                 `B.append` ")" ) pwm
    merge tr = ( B.intercalate "+" $ map _name $ flatten tr
               , dilute $ mergeTreeWeighted align tr)
    tree = buildTree align ms
    align = alignmentBy jsd $ quadPenal gap
{-# INLINE treeMerge #-}

getSuffix :: String -> String
getSuffix = last . splitOn "."

defaultMain :: Options -> IO ()
defaultMain (Options inFl outFl pre th m svg dump gap) = do
    let readMotif = case getSuffix inFl of
                        "fasta" -> readFasta'
                        _ -> readMEME
        writeMotif = case getSuffix outFl of
                        "fasta" -> writeFasta
                        _ -> \fl x -> writeMEME fl x def
    motifs <- readMotif inFl

    if dump
        then
            mapM_ (\(a,b,c) -> B.putStrLn $ B.intercalate "\t" [a, b, B.pack $ show c]) $ pairDistance motifs
        else do
            let motifNumber = length motifs

            hPutStrLn stderr $ printf "Merging Mode: %s" m
            hPutStrLn stderr $ printf "Read %d motifs" motifNumber

            motifs' <- case m of
                "tree" -> do
                    let (newMotifs, tree) = treeMerge th pre motifs gap
                        fn x = B.unpack (_name x) ++ ": " ++ B.unpack (toBS $ toIUPAC $ _pwm x)

                    case svg of
                        Just fl -> do
                            {-
                            let w = 80
                                h = 5 * fromIntegral motifNumber
                            renderCairo fl (dims2D (10*w) (10*h)) $ drawDendrogram w h th tree fn ||| strutX 40
                            -}
                            return newMotifs
                        _ -> return newMotifs
                "iter" -> do
                    let rs = iterativeMerge (alignmentBy jsd (quadPenal gap)) th motifs
                    forM rs $ \(nm, pwm, ws) -> do
                        let pwm' = dilute (pwm, ws)
                        return $ Motif (B.intercalate "+" nm) pwm
                 -- _ -> error "Unkown mode!"

            hPutStrLn stderr $ printf "Write %d motifs" (length motifs')

            writeMotif outFl motifs'
{-# INLINE defaultMain #-}

comb :: [a] -> [(a,a)]
comb (y:ys) = zip (repeat y) ys ++ comb ys
comb _ = []
{-# INLINE comb #-}

main :: IO ()
main = execParser opts >>= defaultMain
  where
    opts = info (helper <*> parser)
            ( fullDesc
           <> header (printf "Merge similar PWMs, version %s" v))
    v = "0.2.0" :: String
