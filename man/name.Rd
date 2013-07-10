\name{name}
\alias{name}
\alias{name<-}
\alias{name,SeqObj-method}
\alias{name<-,SeqObj,character-method}
\alias{name,cobindr-method}
\alias{name<-,cobindr,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  name of cobindR SeqObj object
}
\description{
  name of cobindR seqObj object. 
}
\usage{
\S4method{name}{SeqObj}(x)
\S4method{name}{cobindr}(x)
\S4method{name}{SeqObj,character}(x) <- value
\S4method{name}{cobindr,character}(x) <- value

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR seqObj object
}
  \item{value}{
  the name describing the sequence object
}

}
\value{
  name (character)
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
name(so)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ SeqObj }
\keyword{ cobindR }
