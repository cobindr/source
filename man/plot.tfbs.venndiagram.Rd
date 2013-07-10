\name{plot.tfbs.venndiagram}
\alias{plot.tfbs.venndiagram}
\alias{plot.tfbs.venndiagram-method}
\alias{plot.tfbs.venndiagram,cobindr-method}

\title{
function visualize the overlaps of PWM hits over the sequences.
}
\description{
The distribution of PWM hits over the sequences is visualized as Venn diagram.
If a list of PWM names is provided, only these PWMs are included in the
Venn diagram. If include.empty.seqs == TRUE, sequences without hits of
the specified PWMs are also included in the diagram. If a file path is specified in pdf.name, the diagram will be written into the specified file.
}
\usage{
\S4method{plot.tfbs.venndiagram}{cobindr}(x, pwms, include.empty.seqs = FALSE)
}

\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences and the hits.
}
  \item{pwms}{
    list of PWMs, which shall be visualized in the Venn-Diagram. If no
    list is given, all PWMs in the runObj are used. The package
    "VennDiagram" only allows Venn plots with up to 4 elements. 
}
  \item{include.empty.seqs}{
    logical flag, if include.empty.seqs == TRUE, sequences without hits of
the specified PWMs are also included in the diagram.
}
}

\references{
using the package "VennDiagram" (http://www.biomedcentral.com/1471-2105/12/35/)
}
\author{
Manuela Benary <manuela.benary@cms.hu-berlin.de
}

\section{Warning }{
\itemize{
  \item{unknown identifier:}
  if the list of PWMs contains unknown PWM identifiers a warning is
  given and the method stops
  \item{too many PWMs:}
  if more than 4 PWMs are listed a warning is given and the method stops
  \item{no hits:}
  if no hits are found in the object, the method gives a warning and
  stops
}
} 

\seealso{
\code{\link{plot.tfbs.heatmap}}
}
