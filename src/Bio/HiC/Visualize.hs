{-# LANGUAGE FlexibleContexts #-}
module Bio.HiC.Visualize
    ( drawHiC
    , DrawHiCOpt(..)
    ) where

import qualified Data.Matrix.Generic as MG
import Data.Matrix.Symmetric (SymMatrix(..))
import Codec.Picture
import Data.Colour
import Data.Colour.Names
import Data.Colour.SRGB
import Data.Colour.Palette.BrewerSet
import Data.Default.Class (Default(..))
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U

import Bio.HiC

data DrawHiCOpt = DrawHiCOpt 
    { _range :: !(Maybe (Double, Double))
    , _colours :: !(V.Vector (Colour Double))
    }

instance Default DrawHiCOpt where
    def = DrawHiCOpt
        { _range = Nothing
--        , _colours = V.fromList $ brewerSet Reds 3
        , _colours = V.fromList [white, red]
        }

drawHiC :: FilePath -> ContactMap -> DrawHiCOpt -> IO ()
drawHiC fl m opt = writePng fl $ matToImage (lo,hi) mat $ _colours opt
  where
    mat@(SymMatrix _ v) = _matrix m
    (lo,hi) = case _range opt of
        Just (a,b) -> (a,b)
        _ -> (U.minimum v, U.maximum v)

matToImage :: MG.Matrix m v Double => (Double, Double) -> m v Double -> V.Vector (Colour Double) -> Image PixelRGB8
matToImage (lo,hi) mat cs = generateImage drawPixel r c
  where
    drawPixel i j = pickColour $ mat `MG.unsafeIndex` (i,j)
    (r,c) = MG.dim mat
    pickColour v = colorConvert $ colorMapSmooth (lmap v) cs 
    lmap = linearMapBound (lo,hi) (0,1)
{-# INLINE matToImage #-}

colorConvert :: Colour Double -> PixelRGB8
colorConvert c = let RGB r g b = toSRGB24 c
                 in PixelRGB8 r g b
{-# INLINE colorConvert #-}

-- | map numbers to colors
colorMapSmooth :: Double -- a value from 0 to 1
               -> V.Vector (Colour Double) -> Colour Double
colorMapSmooth x colors = blend p (colors V.! i) $ colors V.! (i+1)
  where
    p = fromIntegral i - x * (fromIntegral n - 1) + 1
    i | x == 1 = n - 2
      | otherwise = truncate $ x * (fromIntegral n - 1)
    n = V.length colors
{-# INLINE colorMapSmooth #-}

linearMapBound :: (Double, Double) -> (Double, Double) -> Double -> Double
linearMapBound (l, u) (l', u') x
    | x < l = l' 
    | x > u = u'
    | otherwise = (x - l) / (u - l) * (u' - l') + l'
{-# INLINE linearMapBound #-}
