\name{write.pairs}
\alias{write.pairs}
\alias{write.pairs-method}
\alias{write.pairs,cobindr-method}

\title{
function to write output of findPairs into file
}

\description{
Function writes the results of findPairs() as a tab-separated file. The file
consists of 6 columns, namely
\itemize{
  \item{a unique id for each pair,}
  \item{the unique id of the sequence, where the pair was found,}
  \item{the names of the corresponding PFMs,}
  \item{the unique id for each PFM, and}
  \item{the distance window in which the pair occurs.}
}}

\usage{
\S4method{write.pairs}{cobindr}(x, file = NULL, background = FALSE)
}
\arguments{
  \item{x}{
    an object of the class "cobindr", which holds all necessary
    information about the sequences and the predicted binding sites.
}
  \item{file}{
    path to file. If filename is 'NULL' a filename is generated based on
    the name of the object of class "cobindr".    
  }
  \item{background}{
    logical flag. If background is 'TRUE' the pairs found in the
    background sequences are used.
  }
}

\author{
Stefan Kroeger <kroeger@informatik.hu-berlin.de>
}

\seealso{
\code{\link{write.bindingsites.table}}, \code{\link{write.bindingsites}}, \code{\link{write.sequences}}, \code{\link{write}}
}

\keyword{file}
\keyword{IO}

