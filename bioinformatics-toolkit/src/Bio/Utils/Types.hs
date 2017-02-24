{-# LANGUAGE GADTs #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE FlexibleContexts #-}
module Bio.Utils.Types
    ( Sorted
    , ordering
    , fromSorted
    , toSorted
    , unsafeToSorted
    ) where

import qualified Data.Foldable as F
import Data.Ord ()

data Sorted f a where
    Sorted :: (F.Foldable f, Ord a)
           => { ordering :: !Ordering
              , fromSorted :: !(f a)
              }
           -> Sorted f a

deriving instance Show (f a) => Show (Sorted f a)

-- | if the data has been sorted, wrap it into Sorted type
toSorted :: (F.Foldable f, Ord a) => f a -> Sorted f a
toSorted xs = Sorted o xs
  where
    o = snd . F.foldl' g (const EQ, EQ) $ xs
    g (func, ord) x
        | ord == EQ = (compare x, ord')
        | ord' == ord || ord' == EQ = (compare x, ord)
        | otherwise = error "data is not sorted"
      where
        ord' = func x

unsafeToSorted :: (F.Foldable f, Ord a) => Ordering -> f a -> Sorted f a
unsafeToSorted = Sorted
