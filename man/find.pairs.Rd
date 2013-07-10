\name{find.pairs}
\alias{find.pairs}

\title{
function to find pairs of binding sites for every sequence in a given
object of class "cobindr"
}
\description{
find.pairs creates a data frame with all pairs in all sequences within
the given distance. 
}
\usage{
find.pairs(x, background_scan = FALSE, n.cpu = NA)
}

\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences and the hits.
}
  \item{background_scan}{
    logical flag, if background_scan = TRUE the pairs for the background
  sequences will be found.
}
\item{n.cpu}{
  number of CPUs to be used for parallelization. Default value is 'NA'
  in which case the number of available CPUs is checked and than used.
  }
}

\value{
  \item{runObj}{an object of the class "cobindr" including the pairs of
  transcription factor binding sites}
}

\author{
Yue-Hien Lee <>
}

\seealso{
  \code{\link{plot.detrending}}
}

\keyword{ manip }
\keyword{ methods }
