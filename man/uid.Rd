\name{uid}
\alias{uid}
\alias{uid<-}
\alias{uid,SeqObj-method}
\alias{uid<-,SeqObj,character-method}
\alias{uid,cobindr-method}
\alias{uid<-,cobindr,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  uid of cobindR SeqObj object
}
\description{
  uid of cobindR seqObj object. 
}
\usage{
\S4method{uid}{SeqObj}(x)
\S4method{uid}{cobindr}(x)
\S4method{uid}{SeqObj,character}(x) <- value
\S4method{uid}{cobindr,character}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR seqObj object
}
  \item{value}{
  the unique id of the sequence or cobindr object
}

}
\value{
  uid (character)
}
\author{
Rob Lehmann  <r.lehmann@biologie.hu-berlin.de>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{uid}},\code{\link{name}},\code{\link{species}},\code{\link{location}},\code{\link{comment}},\code{\link{sequence}}}
\examples{
library(Biostrings)
so <- seqObj(DNAString('A'), id='', name='', species='',comment='',location='')
uid(so)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ SeqObj }
\keyword{ cobindR }
