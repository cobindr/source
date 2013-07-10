\name{plot.pairdistribution}
\alias{plot.pairdistribution}
\alias{plot.pairdistribution-method}
\alias{plot.pairdistribution,cobindr-method}

\title{
function to plot the distribution of the number of pairs in the sequences
}
\description{
For a specified pair of PWMs the function visualizes in how many
sequences how many of the pairs can be found.
}

\usage{
\S4method{plot.pairdistribution}{cobindr}(x, pwm1, pwm2)
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
}

\author{
Manuela Benary <manuela.benary@cms.hu-berlin.de>
}

\seealso{
\code{\link{plot.detrending}}, \code{\link{plot.pairdistance}}
}

\keyword{ dplot }
\keyword{ hplot }
