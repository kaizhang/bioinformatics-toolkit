module Bio.RealWorld.ID where

import qualified Data.ByteString.Char8 as B

class BioID a where
    idToBS :: a -> B.ByteString

newtype UniprotID = UniprotID B.ByteString deriving (Show)

newtype UCSCID = UCSCID B.ByteString deriving (Show)

newtype GOID = GOID B.ByteString deriving (Show)
