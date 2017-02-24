{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE CPP #-}
-- most of the codes in this file are directly copied from JuicyPixel

module Bio.HiC.Visualize.Internal where

#if !MIN_VERSION_base(4,8,0)
import Foreign.ForeignPtr.Safe( ForeignPtr, castForeignPtr )
#else
import Foreign.ForeignPtr( ForeignPtr, castForeignPtr )
#endif

import Foreign.Storable( Storable, sizeOf )
import Control.Monad (when)
import Data.Binary (Binary(..), Get)
import Data.Binary.Get( getWord8
                      , getWord32be
                      , getLazyByteString
                      )
import Data.Binary.Put( runPut
                      , putWord8
                      , putWord32be
                      , putLazyByteString
                      )
import Data.Bits( xor, (.&.), unsafeShiftR )
import qualified Data.Vector.Unboxed as U
import Data.Word
import Data.List (foldl')
import qualified Data.ByteString as B
import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString.Lazy.Char8 as LS
import Data.Vector.Storable ( Vector, unsafeToForeignPtr, unsafeFromForeignPtr0 )
import qualified Data.ByteString.Internal as S
import qualified Data.Vector.Generic as G

import Data.Conduit

-- | Value used to identify a png chunk, must be 4 bytes long.
type ChunkSignature = L.ByteString

-- | Generic header used in PNG images.
data PngIHdr = PngIHdr
    { width             :: !Word32       -- ^ Image width in number of pixel
    , height            :: !Word32       -- ^ Image height in number of pixel
    , bitDepth          :: !Word8        -- ^ Number of bit per sample
    , colourType        :: !PngImageType -- ^ Kind of png image (greyscale, true color, indexed...)
    , compressionMethod :: !Word8        -- ^ Compression method used
    , filterMethod      :: !Word8        -- ^ Must be 0
    , interlaceMethod   :: !PngInterlaceMethod   -- ^ If the image is interlaced (for progressive rendering)
    }
    deriving Show

-- | Data structure during real png loading/parsing
data PngRawChunk = PngRawChunk
    { chunkLength :: Word32
    , chunkType   :: ChunkSignature
    , chunkCRC    :: Word32
    , chunkData   :: L.ByteString
    }

-- | What kind of information is encoded in the IDAT section
-- of the PngFile
data PngImageType =
      PngGreyscale
    | PngTrueColour
    | PngIndexedColor
    | PngGreyscaleWithAlpha
    | PngTrueColourWithAlpha
    deriving Show

-- | Different known interlace methods for PNG image
data PngInterlaceMethod =
      -- | No interlacing, basic data ordering, line by line
      -- from left to right.
      PngNoInterlace

      -- | Use the Adam7 ordering, see `adam7Reordering`
    | PngInterlaceAdam7
    deriving (Enum, Show)

preparePngHeader :: Int -> Int -> PngImageType -> Word8 -> PngIHdr
preparePngHeader w h imgType depth = PngIHdr
    { width             = fromIntegral w
    , height            = fromIntegral h
    , bitDepth          = depth
    , colourType        = imgType
    , compressionMethod = 0
    , filterMethod      = 0
    , interlaceMethod   = PngNoInterlace
    }

instance Binary PngRawChunk where
    put chunk = do
        putWord32be $ chunkLength chunk
        putLazyByteString $ chunkType chunk
        when (chunkLength chunk /= 0)
             (putLazyByteString $ chunkData chunk)
        putWord32be $ chunkCRC chunk

    get = do
        size <- getWord32be
        chunkSig <- getLazyByteString (fromIntegral $ L.length iHDRSignature)
        imgData <- if size == 0
            then return L.empty
            else getLazyByteString (fromIntegral size)
        crc <- getWord32be

        let computedCrc = pngComputeCrc [chunkSig, imgData]
        when (computedCrc `xor` crc /= 0)
             (fail $ "Invalid CRC : " ++ show computedCrc ++ ", "
                                      ++ show crc)
        return PngRawChunk {
            chunkLength = size,
            chunkData = imgData,
            chunkCRC = crc,
            chunkType = chunkSig
        }

instance Binary PngImageType where
    put PngGreyscale = putWord8 0
    put PngTrueColour = putWord8 2
    put PngIndexedColor = putWord8 3
    put PngGreyscaleWithAlpha = putWord8 4
    put PngTrueColourWithAlpha = putWord8 6

    get = get >>= imageTypeOfCode

imageTypeOfCode :: Word8 -> Get PngImageType
imageTypeOfCode 0 = return PngGreyscale
imageTypeOfCode 2 = return PngTrueColour
imageTypeOfCode 3 = return PngIndexedColor
imageTypeOfCode 4 = return PngGreyscaleWithAlpha
imageTypeOfCode 6 = return PngTrueColourWithAlpha
imageTypeOfCode _ = fail "Invalid png color code"

instance Binary PngIHdr where
    put hdr = do
        putWord32be 13
        let inner = runPut $ do
                putLazyByteString iHDRSignature
                putWord32be $ width hdr
                putWord32be $ height hdr
                putWord8    $ bitDepth hdr
                put $ colourType hdr
                put $ compressionMethod hdr
                put $ filterMethod hdr
                put $ interlaceMethod hdr
            crc = pngComputeCrc [inner]
        putLazyByteString inner
        putWord32be crc

    get = do
        _size <- getWord32be
        ihdrSig <- getLazyByteString (L.length iHDRSignature)
        when (ihdrSig /= iHDRSignature)
             (fail "Invalid PNG file, wrong ihdr")
        w <- getWord32be
        h <- getWord32be
        depth <- get
        colorType <- get
        compression <- get
        filtermethod <- get
        interlace <- get
        _crc <- getWord32be
        return PngIHdr {
            width = w,
            height = h,
            bitDepth = depth,
            colourType = colorType,
            compressionMethod = compression,
            filterMethod = filtermethod,
            interlaceMethod = interlace
        }

instance Binary PngInterlaceMethod where
    get = getWord8 >>= \w -> case w of
        0 -> return PngNoInterlace
        1 -> return PngInterlaceAdam7
        _ -> fail "Invalid interlace method"

    put PngNoInterlace    = putWord8 0
    put PngInterlaceAdam7 = putWord8 1

-- signature

-- | Signature signalling that the following data will be a png image
-- in the png bit stream
pngSignature :: ChunkSignature
pngSignature = L.pack [137, 80, 78, 71, 13, 10, 26, 10]

-- | Helper function to help pack signatures.
signature :: String -> ChunkSignature
signature = LS.pack 

-- | Signature for the header chunk of png (must be the first)
iHDRSignature :: ChunkSignature 
iHDRSignature = signature "IHDR"

-- | Signature for a palette chunk in the pgn file. Must
-- occure before iDAT.
pLTESignature :: ChunkSignature
pLTESignature = signature "PLTE"

-- | Signature for a data chuck (with image parts in it)
iDATSignature :: ChunkSignature
iDATSignature = signature "IDAT"

-- | Signature for the last chunk of a png image, telling
-- the end.
iENDSignature :: ChunkSignature
iENDSignature = signature "IEND"

-- | Compute the CRC of a raw buffer, as described in annex D of the PNG
-- specification.
pngComputeCrc :: [L.ByteString] -> Word32
pngComputeCrc = (0xFFFFFFFF `xor`) . L.foldl' updateCrc 0xFFFFFFFF . L.concat
    where updateCrc crc val =
              let u32Val = fromIntegral val
                  lutVal = pngCrcTable U.! fromIntegral ((crc `xor` u32Val) .&. 0xFF)
              in lutVal `xor` (crc `unsafeShiftR` 8)

-- | From the Annex D of the png specification.
pngCrcTable :: U.Vector Word32
pngCrcTable = U.fromListN 256 [ foldl' updateCrcConstant c [zero .. 7] | c <- [0 .. 255] ]
    where zero = 0 :: Int -- To avoid defaulting to Integer
          updateCrcConstant c _ | c .&. 1 /= 0 = magicConstant `xor` (c `unsafeShiftR` 1)
                                | otherwise = c `unsafeShiftR` 1
          magicConstant = 0xedb88320 :: Word32



----------------------------------------------------------------------------

endChunk :: PngRawChunk
endChunk = PngRawChunk { chunkLength = 0
                       , chunkType = iENDSignature
                       , chunkCRC = pngComputeCrc [iENDSignature]
                       , chunkData = L.empty
                       }

type Palette = Vector Word8

preparePalette :: Palette -> PngRawChunk
preparePalette pal = PngRawChunk
  { chunkLength = fromIntegral $ G.length pal
  , chunkType   = pLTESignature
  , chunkCRC    = pngComputeCrc [pLTESignature, binaryData]
  , chunkData   = binaryData
  }
   where binaryData = L.fromChunks [toByteString pal]

toByteString :: forall a. (Storable a) => Vector a -> B.ByteString
toByteString vec = S.PS (castForeignPtr ptr) offset (len * size)
  where (ptr, offset, len) = unsafeToForeignPtr vec
        size = sizeOf (undefined :: a)

generatePng :: Int -> Int -> Palette -> (Int -> Int -> Word8) -> Source m L.ByteString
generatePng w h pal fn = do
    yield pngSignature
    yield $ encode header
    yield $ encode $ preparePalette pal
    loop 0 0
    yield $ encode endChunk
  where
    header = preparePngHeader w h PngIndexedColor 8
    loop i j = do
        fn i j

genericEncodePng :: forall px. (Pixel px, PixelBaseComponent px ~ Word8)
                 => Int -> Int -> Maybe Palette -> PngImageType -> Image px
                 -> L.ByteString
genericEncodePng w h palette imgKind
                 image@(Image { imageWidth = w, imageHeight = h, imageData = arr }) =
  encode PngRawImage { header = hdr
                     , chunks = prependPalette palette [prepareIDatChunk imgEncodedData, endChunk]}
    where hdr = preparePngHeader w h image imgKind 8
          zero = B.singleton 0

          compCount = componentCount (undefined :: px)

          prependPalette Nothing l = l
          prependPalette (Just p) l = preparePalette p : l

          lineSize = compCount * w
          encodeLine line = blitVector arr (line * lineSize) lineSize
          imgEncodedData = Z.compress . L.fromChunks
                        $ concat [[zero, encodeLine line] | line <- [0 .. h - 1]]

prepareIDatChunk :: Lb.ByteString -> PngRawChunk
prepareIDatChunk imgData = PngRawChunk
    { chunkLength = fromIntegral $ Lb.length imgData
    , chunkType   = iDATSignature
    , chunkCRC    = pngComputeCrc [iDATSignature, imgData]
    , chunkData   = imgData
    }

