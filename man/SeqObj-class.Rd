\name{SeqObj-class}
\Rdversion{1.1}
\docType{class}
\alias{SeqObj-class}
\alias{initialize,SeqObj-method}
\alias{rtfbs.intern,SeqObj-method}
\alias{write.fasta,SeqObj-method}

\title{Class \code{"SeqObj"}}
\description{
Container for DNA sequence and its meta-data. 
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SeqObj", seq, id, species, name, comment, location)}.
}

\section{Slots}{
  \describe{
    \item{\code{uid}:}{Object of class \code{"character"} ~~ unique id for
      internal representation}
    \item{\code{name}:}{Object of class \code{"character"} ~~ biological
      reference name, if available}
    \item{\code{species}:}{Object of class \code{"character"} ~~ reference species }
    \item{\code{location}:}{Object of class \code{"character"} ~~ location
      on the reference genome}
    \item{\code{comment}:}{Object of class \code{"character"} ~~ comments
      and notes}
    \item{\code{sequence}:}{Object of class \code{"DNAString"} ~~ the sequence }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "SeqObj")}: ... }
    \item{rtfbs.intern}{\code{signature(object = "SeqObj")}: ... }
    \item{write.fasta}{\code{signature(sequences = "SeqObj")}: ... }
	 }
}

\author{
Manuela Benary <manuela.benary@cms.hu-berlin.de>
}

\seealso{
  \code{\linkS4class{cobindr}}
  \code{\linkS4class{configuration}}
}

\examples{
showClass("SeqObj")
}

\keyword{classes}
