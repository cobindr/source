\name{pValue}
\alias{pValue}
\alias{pValue<-}
\alias{pValue,configuration-method}
\alias{pValue<-,configuration,numeric-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  pValue threshold used for motif hit finding
}
\description{
  pValue threshold used for motif hit finding
}
\usage{
\S4method{pValue}{configuration}(x)
\S4method{pValue}{configuration,numeric}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR configuration object
}
  \item{value}{
  the p-value threshold used for hit searching
}

}
\value{
  pValue threshold (numeric)
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
pValue(cfg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ configuration }
\keyword{ cobindR }
