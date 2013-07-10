\name{configuration}
\alias{configuration}
\alias{configuration<-}
\alias{configuration,cobindr-method}
\alias{configuration<-,cobindr,configuration-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  configuration of cobindr object
}
\description{
  configuration of cobindr object
}
\usage{
\S4method{configuration}{cobindr}(x)
\S4method{configuration}{cobindr,configuration}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindr object
}
  \item{value}{
  returns the configuration object used in this cobindR object
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
\code{\link{uid}},\code{\link{name}},\code{\link{sequences}},\code{\link{bg_sequences}},\code{\link{desc}},\code{\link{configuration}},\code{\link{binding_sites}},\code{\link{bg_binding_sites}},\code{\link{pfm}},\code{\link{pairs}},\code{\link{bg_pairs}},\code{\link{pairs_of_interest}}
}
\examples{

cfg <- cobindRConfiguration()
sequence_type(cfg) <- 'fasta'
sequence_source(cfg) <- system.file('extdata/sox_oct_example_vignette_seqs.fasta',package='cobindR')
sequence_origin(cfg) <- 'Mouse Embryonic Stem Cell Example ChIP-Seq Oct4 Peak configuration'
cbr <- cobindr(cfg)
configuration(cbr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ cobindr }
\keyword{ cobindR }
