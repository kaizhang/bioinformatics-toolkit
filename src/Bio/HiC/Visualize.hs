{-# LANGUAGE FlexibleContexts #-}
module Bio.HiC.Visualize
    ( drawHiC
    , DrawHiCOpt(..)
    ) where

import qualified Data.ByteString.Lazy as L
import qualified Data.Matrix.Generic as MG
import Data.Matrix.Symmetric (SymMatrix(..))
import Codec.Picture
import Data.Colour
import Data.Colour.Names
import Data.Colour.SRGB
import Data.Default.Class (Default(..))
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U

import Bio.HiC

data DrawHiCOpt = DrawHiCOpt 
    { _range :: !(Maybe (Double, Double))
    , _palette :: !(V.Vector (Colour Double))
    }

instance Default DrawHiCOpt where
    def = DrawHiCOpt
        { _range = Nothing
--        , _palette = V.fromList $ brewerSet Reds 3
        , _palette = V.fromList $ interpolate 62 white red
        }

drawHiC :: FilePath -> ContactMap -> DrawHiCOpt -> IO ()
drawHiC fl m opt = case encodePalettedPng pal pic of
    Left e -> error e
    Right png -> L.writeFile fl png
  where
    pic = matToImage nCol (lo,hi) mat
    pal = generateImage fn nCol 1
      where
        fn i _ = colorConvert $ cols V.! i
        cols = _palette opt
    nCol = V.length $ _palette opt
    mat@(SymMatrix _ v) = _matrix m
    (lo,hi) = case _range opt of
        Just (a,b) -> (a,b)
        _ -> (U.minimum v, U.maximum v)

matToImage :: MG.Matrix m v Double => Int -> (Double, Double) -> m v Double -> Image Pixel8
matToImage n (lo,hi) mat = generateImage drawPixel r c
  where
    drawPixel i j | x <= lo = 0
                  | x >= hi = fromIntegral $ n - 1
                  | otherwise = truncate $ (x - lo) / step
      where
        x = mat `MG.unsafeIndex` (i,j)
    step = (hi - lo) / fromIntegral n
    (r,c) = MG.dim mat
{-# INLINE matToImage #-}

colorConvert :: Colour Double -> PixelRGB8
colorConvert c = let RGB r g b = toSRGB24 c
                 in PixelRGB8 r g b
{-# INLINE colorConvert #-}

{-
colorBlend :: Int -> V.Vector (Colour Double) -> V.Vector (Colour Double)
colorBlend n colors | n <= m = colors
                    | otherwise = 
  where
    m = V.length colors
    -}

interpolate :: Int -> Colour Double -> Colour Double -> [Colour Double]
interpolate n c1 c2 = loop 1
  where
    loop i | i > n = []
           | otherwise = blend (fromIntegral i * step) c2 c1 : loop (i+1)
    step = 1 / fromIntegral (n+1)
{-# INLINE interpolate #-}

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
    | isNaN x || x < l = l' 
    | x > u = u'
    | otherwise = (x - l) / (u - l) * (u' - l') + l'
{-# INLINE linearMapBound #-}
