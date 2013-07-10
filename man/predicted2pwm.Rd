\name{predicted2pwm}
\alias{predicted2pwm}
\alias{predicted2pwm-method}
\alias{predicted2pwm,cobindr-method}

\title{
function to convert predicted TFBS hits into a PWM
}
\description{
function converts for each input PWM the predicted TFBS hits into a PWM. Function is intended to be used together with the sequence logo creation function 'plot.tfbslogo'.
}
\usage{
\S4method{predicted2pwm}{cobindr}(x, as.pfm=FALSE)
}

\arguments{
  \item{x}{
    object of class "cobindr" describing the sequences and the
  predicted TFBS.
}
  \item{as.pfm}{logical flag, to indicate whether the function should
  return a PFM (TRUE) or a PWM (FALSE)}
}

\value{
 \item{predPwm}{positional frequency matrix based on consensus matrix}
}


\author{
  Robert Lehmann <r.lehmann@biologie.hu-berlin.de>
}

\seealso{
 \code{\link{plot.tfbslogo}}
}

\keyword{ manip }
\keyword{ methods }

