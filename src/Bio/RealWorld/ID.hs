module Bio.RealWorld.ID where

import qualified Data.ByteString.Char8 as B

class BioID a where
    fromID :: a -> B.ByteString
    toID :: B.ByteString -> a

newtype UniprotID = UniprotID B.ByteString deriving (Show)

newtype UCSCID = UCSCID B.ByteString deriving (Show)

newtype GOID = GOID B.ByteString deriving (Show)

-- | ENCODE Accession
newtype EncodeAcc = EncodeAcc B.ByteString deriving (Show)

-- | Ensembl ID
newtype EnsemblID = EnsemblID B.ByteString deriving (Show)

instance BioID EncodeAcc where
    fromID (EncodeAcc x) = x
    toID = EncodeAcc

instance BioID EnsemblID where
    fromID (EnsemblID x) = x
    toID = EnsemblID
