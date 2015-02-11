{-# LANGUAGE OverloadedStrings #-}
module Bio.GO
    ( GO(..)
    , GOId
    , GOTree
    , getParent
    ) where

import qualified Data.Text as T
import Data.Text.Encoding (decodeUtf8)
import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict as M
import Data.Maybe (fromMaybe)

data GO = GO
    { _oboId :: !GOId
    , _label :: !T.Text
    , _subProcessOf :: !(Maybe GOId)
    , _oboNS :: !T.Text
    }

instance Show GO where
    show (GO a b c d) = T.unpack $ T.intercalate "\t"
        [decodeUtf8 a, b, decodeUtf8 $ fromMaybe "" c, d]

type GOId = B.ByteString

type GOTree = M.HashMap GOId GO

getParent :: GO -> GOTree -> Maybe GO
getParent g tree = case _subProcessOf g of
    Just i -> M.lookup i tree
    _ -> Nothing
{-# INLINE getParent #-}
