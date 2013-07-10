\name{seqObj}
\alias{seqObj}
\alias{seqObj-method}
\alias{seqObj,DNAString,character,character,character,character,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	cobindR SeqObj object constructor
}
\description{
	cobindR SeqObj object constructor
}
\usage{
\S4method{seqObj}{DNAString,character,character,character,character,character}(seq,id,name,species,comment,location)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{
	DNAString object holding the sequence
}
  \item{id}{
	id (character)
}
  \item{name}{
	id (character)
}
  \item{species}{
	id (character)
}
  \item{comment}{
	id (character)
}
  \item{location}{
	id (character)
}
}
\value{
  cobindR SeqObj object
}
\author{
Rob Lehmann  <r.lehmann@biologie.hu-berlin.de>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cobindRConfiguration}}
}
\examples{
	library(Biostrings)
	so <- seqObj(DNAString('A'), id='', name='', species='',comment='',location='')
	sequence(so)

}

\keyword{ constructor}
\keyword{ configuration }
\keyword{ cobindR }
