\name{write.sequences}
\alias{write.sequences}
\alias{write.sequences-method}
\alias{write.sequences,cobindr-method}

\title{
	writes the sequences of a cobindr-object into a fasta file.
}
\description{
	writes the sequences of a cobindr-object into a fasta file.
}
\usage{
\S4method{write.sequences}{cobindr}(x, slotname = "sequences", file = NULL)
}

\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences.
  }
  \item{slotname}{string, describing whether to use foreground
  sequences (default) or background sequences 
  }
  \item{file}{
    path to file. If filename is 'NULL' a filename is generated based on
    the name of the object of class "cobindr".  
}
}

\author{
  Stefan Kroeger <kroeger@informatik.hu-berlin.de>
}

\seealso{
\code{\link{write.bindingsites.table}}, \code{\link{write.bindingsites}},  \code{\link{write.pairs}}, \code{\link{write}}
}

\examples{
cfg <- cobindRConfiguration()
sequence_type(cfg) <- 'fasta'
sequence_source(cfg) <- system.file('extdata/example.fasta', package='cobindR')
# avoid complaint of validation mechanism 
pfm_path(cfg) <- system.file('extdata/pfms',package='cobindR')
pairs(cfg) <- '' 
runObj <- cobindr(cfg)
write.sequences(runObj, file = file.path(tempfile("example.txt",tempdir())) )
}

\keyword{ IO }
\keyword{ file }
