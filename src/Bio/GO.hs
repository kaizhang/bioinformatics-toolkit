{-# LANGUAGE OverloadedStrings #-}
module Bio.GO
    ( GO(..)
    , GOId
    ) where

import qualified Data.Text as T
import Data.Text.Encoding (decodeUtf8)
import qualified Data.ByteString.Char8 as B
import Data.Maybe (fromMaybe)

type GOId = B.ByteString

data GO = GO
    { _label :: !T.Text
    , _oboId :: !GOId
    , _oboNS :: !T.Text
    , _subProcessOf :: !(Maybe GOId)
    }

instance Show GO where
    show (GO a b c d) = T.unpack $ T.intercalate "\t"
        [a, decodeUtf8 b, c, decodeUtf8 $ fromMaybe "" d]
