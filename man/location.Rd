\name{location}
\alias{location}
\alias{location<-}
\alias{location,SeqObj-method}
\alias{location<-,SeqObj,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 location of cobindR SeqObj object
}
\description{
 location of cobindR seqObj object (e.g. chr1)
}
\usage{
\S4method{location}{SeqObj}(x)
\S4method{location}{SeqObj,character}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR seqObj object
}
  \item{value}{
  the location description of the sequence
}

}
\value{
  returns location (character)
}
\author{
Rob Lehmann  <r.lehmann@biologie.hu-berlin.de>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{uid}},\code{\link{name}},\code{\link{species}},\code{\link{location}},\code{\link{comment}},\code{\link{sequence}}
}
\examples{
library(Biostrings)
so <- seqObj(DNAString('A'), id='', name='', species='',comment='',location='')
location(so)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ SeqObj }
\keyword{ cobindR }
