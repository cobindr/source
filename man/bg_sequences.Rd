\name{bg_sequences}
\alias{bg_sequences}
\alias{bg_sequences<-}
\alias{bg_sequences,cobindr-method}
\alias{bg_sequences<-,cobindr,list-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  list of background sequence 
}
\description{
  list of background sequence 
}
\usage{
\S4method{bg_sequences}{cobindr}(x)
\S4method{bg_sequences}{cobindr,list}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindr object
}
  \item{value}{
   list of background sequence of type SeqObj
}
}
\value{
  list of background sequences (SeqObj)
}
\author{
Rob Lehmann  <r.lehmann@biologie.hu-berlin.de>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{uid}},\code{\link{name}},\code{\link{bg_sequences}},\code{\link{bg_sequences}},\code{\link{desc}},\code{\link{configuration}},\code{\link{binding_sites}},\code{\link{bg_binding_sites}},\code{\link{pfm}},\code{\link{pairs}},\code{\link{bg_pairs}},\code{\link{pairs_of_interest}}
}
\examples{

cfg <- cobindRConfiguration()
sequence_type(cfg) <- 'fasta'
sequence_source(cfg) <- system.file('extdata/sox_oct_example_vignette_seqs.fasta',package='cobindR')
sequence_origin(cfg) <- 'Mouse Embryonic Stem Cell Example ChIP-Seq Oct4 Peak bg_sequences'
cbr <- cobindr(cfg)
length(bg_sequences(cbr))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ cobindr }
\keyword{ cobindR }
