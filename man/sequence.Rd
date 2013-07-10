\name{sequence}
\alias{sequence}
\alias{sequence<-}
\alias{sequence,SeqObj-method}
\alias{sequence<-,SeqObj,DNAString-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  returns sequence of cobindR SeqObj object
}
\description{
  returns sequence of cobindR seqObj object. 
}
\usage{
\S4method{sequence}{SeqObj}(x)
\S4method{sequence}{SeqObj,DNAString}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR seqObj object
}
\item{value}{
  DNAString of the actual DNA sequence in this SeqObj
}

}
\value{
  sequence (DNAString)
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
sequence(so)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ SeqObj }
\keyword{ cobindR }
