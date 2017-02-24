{-# LANGUAGE OverloadedStrings #-}
module Bio.GO
    ( GO(..)
    , GOId
    , GOMap
    , getParentById
    , addTerm
    , enrichment
    ) where

import           Bio.Utils.Functions (hyperquick)
import qualified Data.HashMap.Strict as M
import qualified Data.HashSet        as S
import qualified Data.Text           as T

data GO = GO
    { _oboId        :: !GOId
    , _label        :: !T.Text
    , _subProcessOf :: !(Maybe GOId)
    , _oboNS        :: !T.Text
    } deriving (Show, Read)

type GOId = Int

type GOMap = M.HashMap GOId GO

type TermCount = M.HashMap GOId Int

getParentById :: GOId -> GOMap -> Maybe GO
getParentById gid goMap = M.lookup gid goMap >>= _subProcessOf
                                             >>= (`M.lookup` goMap)
{-# INLINE getParentById #-}

-- | Add a GO term to the count table. Term counts will propogate from child to
-- its parents. This function works for cyclical graph as well.
addTerm :: GO -> GOMap -> TermCount -> TermCount
addTerm g m t = loop S.empty g t
  where
    loop visited go table
        | _oboId go `S.member` visited = table
        | otherwise = case _subProcessOf go of
            Nothing -> table'
            Just gid -> loop (S.insert (_oboId go) visited)
                (M.lookupDefault undefined gid m) table'
      where
        table' = M.insertWith (+) (_oboId go) 1 table

enrichment :: (TermCount, Int)  -- ^ Background frequency and the total number
           -> (TermCount, Int)  -- ^ Foreground
           -> [(GOId, Double, Double)]
enrichment (bg, bg_total) (fg, fg_total) =
    flip map (M.toList fg) $ \(gid, fg_count) ->
        let enrich = fromIntegral (fg_count * bg_total) /
                     fromIntegral (fg_total * bg_count)
            bg_count = M.lookupDefault undefined gid bg
            p = 1 - hyperquick fg_count bg_count fg_total bg_total
        in (gid, enrich, p)
