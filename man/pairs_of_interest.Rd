\name{pairs_of_interest}
\alias{pairs_of_interest}
\alias{pairs_of_interest<-}
\alias{pairs_of_interest,cobindr-method}
\alias{pairs_of_interest<-,cobindr,factor-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  pairs_of_interest of cobindr object
}
\description{
  pairs_of_interest of cobindr object. 
}
\usage{
\S4method{pairs_of_interest}{cobindr}(x)
\S4method{pairs_of_interest}{cobindr,factor}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindr object
}
  \item{value}{
  factors specifying the motif pairs that are to be evaluated
}

}
\value{
  pairs_of_interest (factor)
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
sequence_origin(cfg) <- 'Mouse Embryonic Stem Cell Example ChIP-Seq Oct4 Peak pairs_of_interest'
cbr <- cobindr(cfg)
pairs_of_interest(cbr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ cobindr }
\keyword{ cobindR }
