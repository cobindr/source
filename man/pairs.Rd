\name{pairs}
\alias{pairs}
\alias{pairs<-}
\alias{pairs,configuration-method}
\alias{pairs<-,configuration,character-method}
\alias{pairs,cobindr-method}
\alias{pairs<-,cobindr,data.frame-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  motif hit pairs in the foreground sequences
}
\description{
  motif hit pairs in the foreground sequences
}
\usage{
\S4method{pairs}{configuration}(x)
\S4method{pairs}{configuration,character}(x) <- value
\S4method{pairs}{cobindr}(x)
\S4method{pairs}{cobindr,data.frame}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR configuration object
}
  \item{value}{
for a configuration object, pairs is a character specifying the motif pairs which should be considered. for a cobindR object, pairs is a data.frame holding the detected motif pairs.
}

}
\value{
  pairs (character)
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
pairs(cfg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ configuration }
\keyword{ cobindR }
