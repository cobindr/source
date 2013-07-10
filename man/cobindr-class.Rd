\name{cobindr-class}
\Rdversion{1.1}
\docType{class}
\alias{cobindr-class}
\alias{detrending,cobindr-method}
\alias{find.pairs,cobindr-method}
\alias{generate.background,cobindr-method}
\alias{get.bindingsite.ranges,cobindr-method}
\alias{initialize,cobindr-method}
\alias{input.pwm,cobindr-method}
\alias{write,cobindr,character-method}

\title{Class \code{"cobindr"}}
\description{
Container for experiment run and its meta-data
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cobindr", conf, name, desc)}.

}
\section{Slots}{
  \describe{
    \item{\code{uid}:}{Object of class \code{"character"} ~~ unique id for
      internal representation}
    \item{\code{name}:}{Object of class \code{"character"} ~~ name of the experiment}
    \item{\code{sequences}:}{Object of class \code{"list"} ~~ list of
      sequence objects to be analyzed}
    \item{\code{bg_sequences}:}{Object of class \code{"list"} ~~list of
      background sequences for statistical analyses }
    \item{\code{desc}:}{Object of class \code{"character"} ~~ verbal experiment description}
    \item{\code{configuration}:}{Object of class \code{"configuration"} ~~the configuration object used to describe the experiment }
    \item{\code{pfm}:}{Object of class \code{"list"} ~~list of
       pfms to be used  }
    \item{\code{binding_sites}:}{Object of class \code{"data.frame"} ~~ data frame for predicted binding sites. Data frame structure: uid:integer, seqObj_uid:integer, pfm:factor, start:integer, end:integer, score:double, seq:character, strand:factor, source:factor.}
    \item{\code{bg_binding_sites}:}{Object of class \code{"data.frame"} ~~ data frame for predicted binding sites in the background sequences. Data frame structure: uid:integer, seqObj_uid:integer, pfm:factor, start:integer, end:integer, score:double, seq:character, strand:factor, source:factor.}
    \item{\code{pairs}:}{Object of class \code{"data.frame"} ~~ data frame
      for predicted pairs of transcription factors. Data frame
      structure: uid:integer, seqObj_uid:integer, pair:factor,
      bs_uid1:integer, bs_uid2:integer, distance_start:integer.}
    \item{\code{bg_pairs}:}{Object of class \code{"data.frame"} ~~ data frame
      for predicted pairs of transcription factors in the background sequences. Data frame
      structure: uid:integer, seqObj_uid:integer, pair:factor,
      bs_uid1:integer, bs_uid2:integer, distance_start:integer.}
    \item{\code{pairs_of_interest}:}{Object of class \code{"factor"} ~~ contains pairs for search}
  }
}
\section{Methods}{
  \describe{
    \item{detrending}{\code{signature(object = "cobindr")}: ... }
    \item{find.pairs}{\code{signature(object = "cobindr")}: ... }
    \item{generate.background}{\code{signature(object = "cobindr")}: ... }
    \item{get.bindingsite.ranges}{\code{signature(object = "cobindr")}: ... }
    \item{get.pairs}{\code{signature(object = "cobindr")}: ... }
    \item{get.significant.pairs}{\code{signature(object = "cobindr")}: ... }
    \item{initialize}{\code{signature(.Object = "cobindr")}: ... }
    \item{input.pwm}{\code{signature(object = "cobindr")}: ... }
    \item{plot.detrending}{\code{signature(object = "cobindr")}: ... }
    \item{plot.gc}{\code{signature(object = "cobindr")}: ... }
    \item{plot.pairdistance}{\code{signature(object = "cobindr")}: ... }
    \item{plot.pairdistribution}{\code{signature(object = "cobindr")}: ... }
    \item{plot.positionprofile}{\code{signature(object = "cobindr")}: ... }
    \item{plot.positions.simple}{\code{signature(object = "cobindr")}: ... }
    \item{plot.positions}{\code{signature(object = "cobindr")}: ... }
    \item{plot.tfbs.heatmap}{\code{signature(object = "cobindr")}: ... }
    \item{plot.tfbs.venndiagram}{\code{signature(object = "cobindr")}: ... }
    \item{plot.tfbslogo}{\code{signature(object = "cobindr")}: ... }
    \item{predicted2pwm}{\code{signature(object = "cobindr")}: ... }
    \item{rtfbs}{\code{signature(object = "cobindr")}: ... }
    \item{search.gadem}{\code{signature(object = "cobindr")}: ... }
    \item{search.pwm}{\code{signature(object = "cobindr")}: ... }
    \item{testCpG}{\code{signature(object = "cobindr")}: ... }
    \item{write.bindingsites.table}{\code{signature(object = "cobindr")}: ... }
    \item{write.bindingsites}{\code{signature(object = "cobindr")}: ... }
    \item{write.sequences}{\code{signature(object = "cobindr")}: ... }
    \item{write}{\code{signature(x = "cobindr", file = "character")}: ... }
	 }
}

\author{
Manuela Benary <manuela.benary@cms.hu-berlin.de>
}

\seealso{
  \code{\linkS4class{SeqObj}}
  \code{\linkS4class{configuration}}
}

\examples{
showClass("cobindr")
}
\keyword{classes}
