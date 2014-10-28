import qualified Tests.Bed as Bed
import qualified Tests.Motif as Motif
import qualified Tests.ChIPSeq as ChIPSeq
import Test.Tasty

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ Bed.tests
    , ChIPSeq.tests
    , Motif.tests
    ]
