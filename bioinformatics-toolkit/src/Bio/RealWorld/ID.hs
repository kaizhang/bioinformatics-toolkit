module Bio.RealWorld.ID where

import qualified Data.ByteString.Char8 as B

class BioID a where
    fromID :: a -> B.ByteString
    toID :: B.ByteString -> a

newtype UniprotID = UniprotID B.ByteString deriving (Show, Eq)

newtype UCSCID = UCSCID B.ByteString deriving (Show, Eq)

newtype GOID = GOID B.ByteString deriving (Show, Eq)

-- | ENCODE Accession
newtype EncodeAcc = EncodeAcc B.ByteString deriving (Show, Eq)

-- | Ensembl ID
newtype EnsemblID = EnsemblID B.ByteString deriving (Show, Eq)

instance BioID EncodeAcc where
    fromID (EncodeAcc x) = x
    toID = EncodeAcc

instance BioID EnsemblID where
    fromID (EnsemblID x) = x
    toID = EnsemblID
