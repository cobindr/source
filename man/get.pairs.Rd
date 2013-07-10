\name{get.pairs}
\alias{get.pairs}
\alias{get.pairs-method}
\alias{get.pairs,cobindr-method}

\title{
function to get output of findPairs
}

\description{
Function returns the results of findPairs() as a data frame. The data.frame consists of 6 columns, namely
\itemize{
  \item{a unique id for each pair,}
  \item{the unique id of the sequence, where the pair was found,}
  \item{the names of the corresponding PFMs,}
  \item{the unique id for each PFM, and}
  \item{the distance window in which the pair occurs.}
}}

\usage{
\S4method{get.pairs}{cobindr}(x, background = FALSE)
}
\arguments{
  \item{x}{
    an object of the class "cobindr", which holds all necessary
    information about the sequences and the predicted binding sites.
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
\code{\link{get.significant.pairs}}, \code{\link{write.bindingsites}}, \code{\link{write.sequences}}, \code{\link{write}}
}

\keyword{IO}

