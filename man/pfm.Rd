\name{pfm}
\alias{pfm}
\alias{pfm<-}
\alias{pfm,cobindr-method}
\alias{pfm<-,cobindr,list-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  pfm list used in experiment
}
\description{
  pfm list used in experiment
}
\usage{
\S4method{pfm}{cobindr}(x)
\S4method{pfm}{cobindr,list}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindr object
}
  \item{value}{
  a list of motif matrices
}

}
\value{
  pfm (list of motif matrices)
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
sequence_origin(cfg) <- 'Mouse Embryonic Stem Cell Example ChIP-Seq Oct4 Peak pfm'
cbr <- cobindr(cfg)
pfm(cbr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ cobindr }
\keyword{ cobindR }
