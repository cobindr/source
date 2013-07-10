\name{plot.tfbs.heatmap}
\alias{plot.tfbs.heatmap}
\alias{plot.tfbs.heatmap-method}
\alias{plot.tfbs.heatmap,cobindr-method}

\title{
function to do plot a heatmap of overlaps between all specified PWMs
}
\description{
plot.tfbs.heatmap plots a heatmap of overlaps between all specified PWMs. For
each overlap, the significance is determined based on the hypergeometric
test. If a file path is specified in pdf.name, the diagram will be written into the specified file.
}
\usage{
\S4method{plot.tfbs.heatmap}{cobindr}(x, pwms, include.empty.seqs = FALSE)
}

\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences and the hits.
}
  \item{pwms}{
    list of PWMs, for which the overlap will be visualized. If no list
  is given, all PWMs in runObj are used.
}
  \item{include.empty.seqs}{
    logical flag, if include.empty.seqs == TRUE, sequences without hits of
the specified PWMs are also included in the diagram.
}
}

\details{
  In this plot for each pair of PWMs the overlap of sequences with hits
  of the given PWMs is calculated. The number of sequences in each
  overlap are color-coded in the heatmap. For each overlap the
  significance is calculated using the hypergeometric test. If the
  significance is below 0.05 (or below 0.01), the corresponding field is
  marked with one (or two) *.
}

\author{
Manuela Benary <manuela.benary@cms.hu-berlin.de>
}

\section{Warning }{
\itemize{
  \item{unknown identifier}
  if the list of PWMs contains unknown PWM identifiers a warning is
  given and the method stops
  \item{no hits}
  if no hits are found in the object, the method gives a warning and
  stops
}
} 

\seealso{
\code{\link{plot.tfbs.venndiagram}}
}

