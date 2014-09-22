import Bio.Motif
import Bio.Utils.Fasta
import Data.Conduit
import qualified Data.Conduit.List as CL

main = do x <- readFasta "motifs.fasta" $$ CL.consume :: IO [Motif]
          print x
