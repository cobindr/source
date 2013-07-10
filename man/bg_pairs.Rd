\name{bg_pairs}
\alias{bg_pairs}
\alias{bg_pairs<-}
\alias{bg_pairs,cobindr-method}
\alias{bg_pairs<-,cobindr,data.frame-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  motif hit pairs in the background sequences
}
\description{
  motif hit pairs in the background sequences
}
\usage{
\S4method{bg_pairs}{cobindr}(x)
\S4method{bg_pairs}{cobindr,data.frame}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindr object
}
  \item{value}{
  data.frame holding the binding site pairs in the background sequences
}

}
\value{
  background motif pairs (data.frame)
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
sequence_origin(cfg) <- 'Mouse Embryonic Stem Cell Example ChIP-Seq Oct4 Peak bg_pairs'
cbr <- cobindr(cfg)
bg_pairs(cbr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ cobindr }
\keyword{ cobindR }
