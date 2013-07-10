\name{sequences}
\alias{sequences}
\alias{sequences<-}
\alias{sequences,cobindr-method}
\alias{sequences<-,cobindr,list-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  sequences of cobindr object
}
\description{
  sequences of cobindr object. 
}
\usage{
\S4method{sequences}{cobindr}(x) 
\S4method{sequences}{cobindr,list}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindr object
}
  \item{value}{
  the list of input sequences of type SeqObj
}

}
\value{
  sequences (character)
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
sequence_origin(cfg) <- 'Mouse Embryonic Stem Cell Example ChIP-Seq Oct4 Peak Sequences'
cbr <- cobindr(cfg)
length(sequences(cbr))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ cobindr }
\keyword{ cobindR }
