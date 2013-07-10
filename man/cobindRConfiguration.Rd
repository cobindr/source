\name{cobindRConfiguration}
\alias{cobindRConfiguration}
\alias{cobindRConfiguration-method}
\alias{cobindRConfiguration,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	cobindR configuration object constructor
}
\description{
	cobindR configuration object constructor
}
\usage{
\S4method{cobindRConfiguration}{character}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	path to configuration file. NULL by default
}
}
\value{
  cobindR configuration object
}
\author{
Rob Lehmann  <r.lehmann@biologie.hu-berlin.de>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{seqObj}}
}
\examples{
cfg <- cobindRConfiguration()
}

\keyword{ constructor}
\keyword{ configuration }
\keyword{ cobindR }
