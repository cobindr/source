\name{plot.positions.simple}
\alias{plot.positions.simple}
\alias{plot.positions.simple-method}
\alias{plot.positions.simple,cobindr-method}

\title{
function to plot  hits for each PWM on the individual sequence
}
\description{
plot.positions plots hits for each PWM on the individual sequence. Which
sequences to plot can be specified by providing a list of sequence
identifiers seq.ids. Which PWMs to plot can be specified as list of
PWMs. The total height of the plot can be adjusted via argument height.
}
\usage{
\S4method{plot.positions.simple}{cobindr}(x, seq.ids, pwms, main)
}

\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences and the hits.
}
  \item{seq.ids}{
    list of sequence identifiers, for which the positions of TFBS will
    be plotted.
}
  \item{pwms}{
    list of PWMs, for which the positions will be visualized. If no list
    is given, all PWMs in runObj are used.
}
  \item{main}{
    title for the plot, if no title is given than 'predicted TFBS
    positions per sequence' will be used
}

}

\author{
Robert Lehmann <r.lehmann@biologie.hu-berlin.de
}


\seealso{
\code{\link{plot.positionprofile}}
}

\keyword{ dplot }
\keyword{ hplot }
