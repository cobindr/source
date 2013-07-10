\name{fdrThreshold}
\alias{fdrThreshold}
\alias{fdrThreshold<-}
\alias{fdrThreshold,configuration-method}
\alias{fdrThreshold<-,configuration,numeric-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  fdrThreshold of cobindR configuration object
}
\description{
  fdrThreshold of cobindR configuration object. 
}
\usage{
\S4method{fdrThreshold}{configuration}(x)
\S4method{fdrThreshold}{configuration,numeric}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR configuration object
}
  \item{value}{
  the false discovery rate threshold to be used for hit search
}

}
\value{
  fdrThreshold (numeric)
}
\author{
Rob Lehmann  <r.lehmann@biologie.hu-berlin.de>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{id}},\code{\link{experiment_description}},\code{\link{sequence_source}},\code{\link{sequence_origin}},\code{\link{sequence_type}},\code{\link{bg_sequence_source}},\code{\link{bg_sequence_origin}},\code{\link{bg_sequence_type}},\code{\link{species}},\code{\link{downstream}},\code{\link{upstream}},\code{\link{max_distance}},\code{\link{pairs}},\code{\link{pfm_path}},\code{\link{threshold}},\code{\link{fdrThreshold}},\code{\link{path}},\code{\link{mart}},\code{\link{pValue}}
}
\examples{
cfg <- cobindRConfiguration()
fdrThreshold(cfg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ configuration }
\keyword{ cobindR }
