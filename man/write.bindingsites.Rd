\name{write.bindingsites}
\alias{write.bindingsites}
\alias{write.bindingsites-method}
\alias{write.bindingsites,cobindr-method}

\title{
	writes predicted binding sites as a BED file.
}
\description{
	writes predicted binding sites as a BED file.
}
\usage{
\S4method{write.bindingsites}{cobindr}(x, file = NULL, background = FALSE)
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
    logical flag. If background is 'TRUE' the binding sites found in the
    background sequences are used.
  }
}

\author{
Stefan Kroeger <kroeger@informatik.hu-berlin.de>
}
\note{
At the moment write.bindingsites() only works for sequences based on gene
ids. Otherwise please use write.bindingsites.table().
}

\seealso{
\code{\link{write.bindingsites.table}}, \code{\link{write.pairs}}, \code{\link{write.sequences}}, \code{\link{write}}
}

\keyword{file}
\keyword{IO}