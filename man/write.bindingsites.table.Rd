\name{write.bindingsites.table}
\alias{write.bindingsites.table}
\alias{write.bindingsites.table-method}
\alias{write.bindingsites.table,cobindr-method}

\title{
function to write predicted TFBS into a tab-separated file.
}
\description{
function to write predicted TFBS into a tab-separated file.
}
\usage{
\S4method{write.bindingsites.table}{cobindr}(x, file = NULL)
}

\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences and the predicted binding sites.
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
\code{\link{write.pairs}}, \code{\link{write.bindingsites}}, \code{\link{write.sequences}}, \code{\link{write}}
}

\keyword{file}
\keyword{IO}
