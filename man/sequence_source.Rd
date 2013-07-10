\name{sequence_source}
\alias{sequence_source}
\alias{sequence_source<-}
\alias{sequence_source,configuration-method}
\alias{sequence_source<-,configuration,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  returns sequence_source of cobindR configuration object
}
\description{
  returns sequence_source of cobindR configuration object. 
}
\usage{
\S4method{sequence_source}{configuration}(x)
\S4method{sequence_source}{configuration,character}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR configuration object
}
  \item{value}{
  the source of which the sequence is retrieved
}

}
\value{
  sequence_source (character)
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
sequence_source(cfg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ configuration }
\keyword{ cobindR }
