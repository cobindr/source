\name{plot.positions}
\alias{plot.positions}
\alias{plot.positions-method}
\alias{plot.positions,cobindr-method}

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
\S4method{plot.positions}{cobindr}(x, seq.ids, pwms, main, order.seq = FALSE, wind.size = 400, frac = 10)
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
  \item{order.seq}{
    logical flag, if TRUE similar patterns of TFBS are shown
    together. This is computationally expensive for large numbers of sequences.
}
  \item{wind.size}{
    integer describing the windows which will be used to enhance
    clustering of TFBS patterns. Necessary if order.seq=TRUE
}
  \item{frac}{
    integer
}
}

\author{
Robert Lehmann <r.lehmann@biologie.hu-berlin.de
}


%\seealso{
%\code{\link{plot.positionprofile}}
%}

\keyword{ dplot }
\keyword{ hplot }

