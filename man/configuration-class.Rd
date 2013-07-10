\name{configuration-class}
\Rdversion{1.1}
\docType{class}
\alias{configuration-class}
\alias{initialize,configuration-method}
\alias{read.background.fasta,configuration-method}
\alias{read.pfm,configuration-method}
\alias{read.sequences,configuration-method}
\alias{write,configuration,character-method}

\title{Class \code{"configuration"}}
\description{
Container for experiment description.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("configuration", fname)}.
}
\section{Slots}{
  \describe{
    \item{\code{id}:}{Object of class \code{"character"} ~~ unique id for
      internal representation }
    \item{\code{experiment_description}:}{Object of class \code{"character"} ~~ verbal experiment description}
    \item{\code{sequence_source}:}{Object of class \code{"character"} ~~ file path or list of paths}
    \item{\code{sequence_origin}:}{Object of class \code{"character"} ~~ source of sequence data, e.g. ensembl}
    \item{\code{sequence_type}:}{Object of class \code{"character"} ~~ either ChipSeq or Fasta or BED are available}
    \item{\code{bg_sequence_source}:}{Object of class \code{"character"} ~~ file path or list of paths}
    \item{\code{bg_sequence_origin}:}{Object of class \code{"character"}
      ~~ how the background is obtained - either simulated or from fasta files or from gene ids }
    \item{\code{bg_sequence_type}:}{Object of class \code{"character"}
      ~~ determines the generation of the background sequences. Possible values: simulated, fasta and geneid}
    \item{\code{species}:}{Object of class \code{"character"} ~~ reference species }
    \item{\code{downstream}:}{Object of class \code{"numeric"} ~~ length of
    sequence downstream of reference point, e.g. transcription start site}
    \item{\code{upstream}:}{Object of class \code{"numeric"} ~~ length of
      sequence upstream of reference point, e.g. transcription start site }
    \item{\code{max_distance}:}{Object of class \code{"numeric"} ~~ maximal
      distance allowed between cooccuring transcription factor binding sites}
    \item{\code{pairs}:}{Object of class \code{"character"} ~~ list of
      pairs of interesting transcription factors }
    \item{\code{pfm_path}:}{Object of class \code{"character"} ~~ path
      to pfm matrix file}
    \item{\code{threshold}:}{Object of class \code{"numeric"} ~~ threshold
      for transcription factor binding site prediction}
    \item{\code{fdrThreshold}:}{Object of class \code{"numeric"} ~~ false
      discovery rate for filtering results (used in rtfbs)}
    \item{\code{date}:}{Object of class \code{"character"} ~~ data of
      experiment run }
    \item{\code{path}:}{Object of class \code{"character"} ~~ path of
      configuration file}
    \item{\code{mart}:}{Object of class \code{"character"} ~~ optional
      mirror for biomart}
    \item{\code{pseudocount}:}{Object of class \code{"numeric"} ~~ sets the pseudocount for
      the detrending analysis}
    \item{\code{pValue}:}{Object of class \code{"numeric"} ~~ optional
      p-Value for search with RGadem}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "configuration")}: ... }
    \item{read.background.fasta}{\code{signature(object = "configuration")}: ... }
    \item{read.pfm}{\code{signature(object = "configuration")}: ... }
    \item{read.sequences}{\code{signature(object = "configuration")}: ... }
    \item{write}{\code{signature(x = "configuration", file = "character")}: ... }
	 }
}

\author{
Manuela Benary <manuela.benary@cms.hu-berlin.de>
}

\seealso{
  \code{\linkS4class{SeqObj}}
  \code{\linkS4class{cobindr}}
}

\examples{
showClass("configuration")
}
\keyword{classes}
