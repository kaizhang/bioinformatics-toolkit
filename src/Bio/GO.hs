{-# LANGUAGE OverloadedStrings #-}
module Bio.GO
    ( GO(..)
    , GOId
    , GOMap
    , getParentById
    ) where

import qualified Data.Text as T
import Data.Text.Encoding (decodeUtf8)
import qualified Data.ByteString.Char8 as B
import qualified Data.HashMap.Strict as M
import Data.Maybe

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

type GOMap = M.HashMap GOId GO

getParentById :: GOId -> GOMap -> Maybe GO
getParentById gid goMap = M.lookup gid goMap >>= _subProcessOf
                                             >>= (`M.lookup` goMap)
{-# INLINE getParentById #-}

{-
buildGOTree :: GOMap -> [Tree GO]
buildGOTree goMap = map build roots
  where
    build = unfoldTree $ \x -> (x, M.lookupDefault [] (_oboId x) childrenMap)

    childrenMap = foldl' f M.empty goMap
      where
        f m go = case _subProcessOf go of
            Just x -> M.insertWith (++) x [go] m
            _ -> m
    roots = mapMaybe f . M.keys $ childrenMap
      where
        f x = do g <- M.lookup x goMap
                 if isNothing $ _subProcessOf g
                    then Nothing
                    else Just g
                    -}
