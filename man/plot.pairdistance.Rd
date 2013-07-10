\name{plot.pairdistance}
\alias{plot.pairdistance}
\alias{plot.pairdistance-method}
\alias{plot.pairdistance,cobindr-method}

\title{
function to plot the distance of the pairs in the sequences
}
\description{
For a specified pair of PWMs the function creates histogram plot of distances between pairs of TFs as specified by pwm1 and pwm2
}

\usage{
\S4method{plot.pairdistance}{cobindr}(x, pwm1, pwm2, breaks=50, main=NA, xlab=NA, ylab=NA, background=FALSE)
}
\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences and the hits.
}
  \item{pwm1}{
    name of the first PWM

}
  \item{pwm2}{
    name of the second PWM
}
  \item{breaks}{ number of breaks to separate the distance distribution into}
  \item{main}{figure title}
  \item{xlab}{label for the x-axis of the figure}
  \item{ylab}{label for the y-axis of the figure}
  \item{background}{flag allowing to plot foreground or background distance distribution}
}

\author{
Manuela Benary <manuela.benary@cms.hu-berlin.de>
}

\seealso{
\code{\link{plot.pairdistribution}}
}

\keyword{ dplot }
\keyword{ hplot }
