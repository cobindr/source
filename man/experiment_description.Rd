\name{experiment_description}
\alias{experiment_description}
\alias{experiment_description<-}
\alias{experiment_description,configuration-method}
\alias{experiment_description<-,configuration,character-method}
\alias{experiment_description,cobindr-method}
\alias{experiment_description<-,cobindr,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  description of cobindR or configuration object
}
\description{
  description of cobindR or configuration object
}
\usage{
\S4method{experiment_description}{configuration}(x)
\S4method{experiment_description}{configuration,character}(x) <- value
\S4method{experiment_description}{cobindr}(x)
\S4method{experiment_description}{cobindr,character}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR or configuration object
}
  \item{value}{
  	description 
}

}
\value{
  experiment description (character)
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

experiment_description(cfg)

sequence_type(cfg) <- 'fasta'
sequence_source(cfg) <- system.file('extdata/sox_oct_example_vignette_seqs.fasta',package='cobindR')
sequence_origin(cfg) <- 'Mouse Embryonic Stem Cell Example ChIP-Seq Oct4 Peak desc'
cbr <- cobindr(cfg)

experiment_description(cbr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ configuration }
\keyword{ cobindR }
