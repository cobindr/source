\name{get.significant.pairs}
\alias{get.significant.pairs}
\alias{get.significant.pairs-method}
\alias{get.significant.pairs,cobindr-method}

\title{
function to returns the results of detrending as a data.frame
}
\description{
get.significant.pairs returns a data.frame of observed distances between the specified pair of PWMs in the foreground set of the sequences as well as the background set of sequences. The distance distribution for the pair
in the background is used for detrending. 
}

\usage{
\S4method{get.significant.pairs}{cobindr}(x, pwm1, pwm2, bin_length=20, z_value=3, overlap=0, abs.distance=FALSE)
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
  \item{bin_length}{
    defines size of bins for distance analysis, default value is 20nucleotides
}
\item{z_value}{level of significance}
\item{overlap}{number of nucleotides which are allowed for an overlap}
\item{abs.distance}{logical flag}
}

\author{
Stefan Kroeger <kroeger@informatik.hu-berlin.de>
}

\seealso{
\code{\link{plot.detrending}}, \code{\link{get.pairs}}, \code{\link{find.pairs}}
}

\keyword{IO}
\keyword{pairs}

