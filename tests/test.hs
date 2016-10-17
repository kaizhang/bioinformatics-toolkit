import qualified Tests.Bed as Bed
import qualified Tests.Motif as Motif
import qualified Tests.Seq as Seq
import qualified Tests.GREAT as GREAT
import qualified Tests.Tools as Tools
import Test.Tasty

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ Bed.tests
    , Seq.tests
    , Motif.tests
    , GREAT.tests
    , Tools.tests
    ]
