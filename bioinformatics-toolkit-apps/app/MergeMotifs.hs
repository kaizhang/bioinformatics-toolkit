{-# LANGUAGE OverloadedStrings #-}

module Main where

import           AI.Clustering.Hierarchical        hiding (drawDendrogram)
import           Control.Monad                     (forM, forM_)
import qualified Data.ByteString.Char8             as B
import           Data.Default
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.List.Split                   (splitOn)
import           Data.Monoid                       ((<>))
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
import           Bio.Seq                           (toBS)
import           Bio.Utils.Functions

data CMD = Merge { mergeInput  :: FilePath
                 , mode        :: String
                 , thres       :: Double
                 , alignOption :: AlignOption
                 , mergeOutput :: FilePath
                 }
         | Dist { inputA      :: FilePath
                , inputB      :: Maybe FilePath
                , alignOption :: AlignOption
                }
         | Cluster { clusterInput :: !FilePath
                   , cutoff       :: !Double
                   , alignOption  :: !AlignOption
                   }

mergeParser :: Parser CMD
mergeParser = Merge
    <$> argument str (metavar "INPUT")
    <*> strOption
        ( long "mode"
       <> short 'm'
       <> value "iter"
       <> metavar "MODE"
       <> help "Merging algorithm, could be iter or tree, default is iter" )
    <*> option auto
        ( long "thres"
       <> short 't'
       <> value 0.2
       <> metavar "THRESHOLD"
       <> help "two motifs that have distance belowing threshold would be merged, default is 0.2" )
    <*> alignParser
    <*> strOption
        ( long "output"
       <> short 'o'
       <> value "merged_output.meme"
       <> metavar "OUTPUT" )

distParser :: Parser CMD
distParser = Dist
    <$> strOption
           ( short 'a'
          <> metavar "INPUT_A" )
    <*> (optional . strOption)
           ( short 'b'
          <> metavar "INPUT_B" )
    <*> alignParser

clusterParser :: Parser CMD
clusterParser = Cluster
    <$> argument str (metavar "INPUT")
    <*> option auto
        ( long "height"
       <> short 'h'
       <> value 0.2
       <> metavar "HEIGHT"
       <> help "Cut hierarchical tree at given height. Default: 0.2" )
    <*> alignParser

data AlignOption = AlignOption
    { gap     :: Double
    , gapMode :: String
    , avgMode :: String
    }

alignParser :: Parser AlignOption
alignParser = AlignOption
     <$> option auto
           ( long "gap"
          <> short 'g'
          <> value 0.05
          <> metavar "GAP_PENALTY"
          <> help "Gap penalty, default: 0.05" )
     <*> strOption
           ( long "gap_mode"
          <> value "exp"
          <> metavar "GAP_MODE"
          <> help "Gap penalty mode, one of linear, quadratic, cubic, and exp. default: exp." )
     <*> strOption
           ( long "avg_mode"
          <> value "l1"
          <> metavar "AVERAGE_MODE"
          <> help "Averaging function, one of l1, l2, l3, max. default: l1." )

treeMerge :: Double -> String -> [Motif] -> AlignOption
          -> ([Motif], Dendrogram Motif)
treeMerge th pre ms alignOpt = (zipWith f [0::Int ..] $ map merge $ tree `cutAt` th, tree)
  where
    f i (suffix, pwm) = Motif ((B.pack $ pre ++ "_" ++ show i ++ "_" ++ show (toIUPAC pwm))
                                 `B.append` "("
                                 `B.append` suffix
                                 `B.append` ")" ) pwm
    merge tr = ( B.intercalate "+" $ map _name $ flatten tr
               , dilute $ mergeTreeWeighted align tr)
    tree = buildTree align ms
    align = mkAlignFn alignOpt
{-# INLINE treeMerge #-}

getSuffix :: String -> String
getSuffix = last . splitOn "."
{-# INLINE getSuffix #-}

mkAlignFn :: AlignOption -> AlignFn
mkAlignFn (AlignOption gap gapMode avgMode) = alignmentBy jsd (pFn gap) avgFn
  where
    pFn = case gapMode of
        "linear"    -> linPenal
        "quadratic" -> quadPenal
        "cubic"     -> cubPenal
        "exp"       -> expPenal
        _           -> error "Unknown gap mode"
    avgFn = case avgMode of
        "l1"  -> l1
        "l2"  -> l2
        "l3"  -> l3
        "max" -> lInf
        _     -> error "Unknown average mode"
{-# INLINE mkAlignFn #-}

readMotif :: FilePath -> IO [Motif]
readMotif fl = case getSuffix fl of
    "fasta" -> readFasta' fl
    _       -> readMEME fl

writeMotif :: FilePath -> [Motif] -> IO ()
writeMotif fl motifs = case getSuffix fl of
    "fasta" -> writeFasta fl motifs
    _       -> writeMEME fl motifs def


defaultMain :: CMD -> IO ()
defaultMain (Dist a b alignOpt) = do
    motifsA <- readMotif a
    pairs <- case b of
        Nothing -> return $ comb motifsA
        Just b' -> do
            motifsB <- readMotif b'
            return [(x,y) | x <- motifsA, y <- motifsB]
    forM_ pairs $ \(x,y) -> do
        let (d,_) = alignFn (_pwm x) $ (_pwm y)
        B.putStrLn $ B.intercalate "\t" [_name x, _name y, toShortest d]
  where
    alignFn = mkAlignFn alignOpt

defaultMain (Merge inFl m th alignOpt outFl)= do
    motifs <- readMotif inFl
    let motifNumber = length motifs

    hPutStrLn stderr $ printf "Merging Mode: %s" m
    hPutStrLn stderr $ printf "Read %d motifs" motifNumber

    motifs' <- case m of
        "tree" -> do
            let (newMotifs, tree) = treeMerge th "" motifs alignOpt
                fn x = B.unpack (_name x) ++ ": " ++ B.unpack (toBS $ toIUPAC $ _pwm x)

            {-
            case svg of
                Just fl -> do
                    let w = 80
                        h = 5 * fromIntegral motifNumber
                    renderCairo fl (dims2D (10*w) (10*h)) $ drawDendrogram w h th tree fn ||| strutX 40
                    return newMotifs
                    -}
            return newMotifs
        "iter" -> do
            let rs = iterativeMerge (mkAlignFn alignOpt) th motifs
            forM rs $ \(nm, pwm, ws) -> do
                let pwm' = dilute (pwm, ws)
                return $ Motif (B.intercalate "+" nm) pwm
         -- _ -> error "Unkown mode!"

    hPutStrLn stderr $ printf "Write %d motifs" (length motifs')

    writeMotif outFl motifs'

defaultMain (Cluster inFl c alignOpt) = do
    motifs <- readMotif inFl
    let tree = buildTree align motifs
        align = mkAlignFn alignOpt
    forM_ (tree `cutAt` c) $ \t ->
        B.putStrLn $ B.intercalate "\t" $ map _name $ flatten t
{-# INLINE defaultMain #-}

comb :: [a] -> [(a,a)]
comb (y:ys) = zip (repeat y) ys ++ comb ys
comb _      = []
{-# INLINE comb #-}

main :: IO ()
main = execParser opts >>= defaultMain
  where
    opts = info (helper <*> parser)
            ( fullDesc
           <> header (printf "Compare, align and merge motifs, version %s" v))
    v = "0.2.0" :: String
    parser = subparser $ (
        command "merge" (info (helper <*> mergeParser) $
            fullDesc <> progDesc "Merge motifs")
     <> command "dist" (info (helper <*> distParser) $
            fullDesc <> progDesc "Align and compute pairwise distances")
     <> command "cluster" (info (helper <*> clusterParser) $
            fullDesc <> progDesc "Perform hierarchical clustering on motifs")
     )
