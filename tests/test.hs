import qualified Tests.Bed as Bed
import qualified Tests.Bam as Bam
import qualified Tests.Motif as Motif
import qualified Tests.ChIPSeq as ChIPSeq
import qualified Tests.Seq as Seq
import qualified Tests.GREAT as GREAT
import Test.Tasty

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ Bed.tests
    , Bam.tests
    , Seq.tests
    , ChIPSeq.tests
    , Motif.tests
    , GREAT.tests
    ]
