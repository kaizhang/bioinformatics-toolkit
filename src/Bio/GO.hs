module Bio.GO
    ( GO(..)
    , GOId
    ) where

import qualified Data.Text as T
import qualified Data.ByteString.Char8 as B

type GOId = B.ByteString

data GO = GO
    { _label :: !T.Text
    , _oboId :: !GOId
    , _oboNS :: !T.Text
    } deriving (Show)
