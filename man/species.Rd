\name{species}
\alias{species}
\alias{species<-}
\alias{species,configuration-method}
\alias{species<-,configuration,character-method}
\alias{species,SeqObj-method}
\alias{species<-,SeqObj,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  species of cobindR configuration or SeqObj
}
\description{
  species of cobindR configuration or SeqObj
}
\usage{
\S4method{species}{configuration}(x)
\S4method{species}{configuration,character}(x) <- value
\S4method{species}{SeqObj}(x)
\S4method{species}{SeqObj,character}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR configuration object
}
  \item{value}{
  name of species in this experiment or SeqObj
}

}
\value{
  sequence / experiment species (character)
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
species(cfg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ configuration }
\keyword{ cobindR }
