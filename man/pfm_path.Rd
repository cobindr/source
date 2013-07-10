\name{pfm_path}
\alias{pfm_path}
\alias{pfm_path<-}
\alias{pfm_path,configuration-method}
\alias{pfm_path<-,configuration,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  path to pfms to be used 
}
\description{
  path to pfms to be used 
}
\usage{
\S4method{pfm_path}{configuration}(x)
\S4method{pfm_path}{configuration,character}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR configuration object
}
  \item{value}{
  the path to the folder containing the motif matrices to be used 
}

}
\value{
  pfm_path (character)
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
pfm_path(cfg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ configuration }
\keyword{ cobindR }
