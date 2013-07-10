\name{comment}
\alias{comment}
\alias{comment<-}
\alias{comment,SeqObj-method}
\alias{comment<-,SeqObj,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  comment of cobindR SeqObj object
}
\description{
  comment of cobindR SeqObj object
}
\usage{
\S4method{comment}{SeqObj}(x)
\S4method{comment}{SeqObj,character}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a cobindR seqObj object
	}
	  \item{value}{
	  comment to the sequence (character)
	}
}

\value{
  comment (character)
}
\author{
Rob Lehmann  <r.lehmann@biologie.hu-berlin.de>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{uid}},\code{\link{name}},\code{\link{species}},\code{\link{comment}},\code{\link{location}},\code{\link{sequence}}
}
\examples{
library(Biostrings)
	so <- seqObj(DNAString('A'), id='', name='', species='',comment='',location='')
comment(so)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ accessor}
\keyword{ SeqObj }
\keyword{ cobindR }
