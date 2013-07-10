\name{pseudocount}
\alias{pseudocount}
\alias{pseudocount<-}
\alias{pseudocount,configuration-method}
\alias{pseudocount<-,configuration,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  pseudocount of cobindR configuration object
}
\description{
  pseudocount of cobindR configuration object. Set to 10 as default
}
\usage{
\S4method{pseudocount}{configuration}(x)
\S4method{pseudocount}{configuration,character}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR configuration object
}
  \item{value}{
  pseudocount for detrending analysis, i.e. the default number in each distance bin.
}

}
\value{
  pseudocount (numeric)
}
\author{
Rob Lehmann  <r.lehmann@biologie.hu-berlin.de>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{id}},\code{\link{experiment_description}},\code{\link{sequence_source}},\code{\link{sequence_origin}},\code{\link{sequence_type}},\code{\link{bg_sequence_source}},\code{\link{bg_sequence_origin}},\code{\link{bg_sequence_type}},\code{\link{species}},\code{\link{downstream}},\code{\link{upstream}},\code{\link{max_distance}},\code{\link{pairs}},\code{\link{pfm_path}},\code{\link{threshold}},\code{\link{fdrThreshold}},\code{\link{path}},\code{\link{mart}},\code{\link{pseudocount}},\code{\link{pValue}}
}
\examples{
cfg <- cobindRConfiguration()
pseudocount(cfg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ configuration }
\keyword{ cobindR }
