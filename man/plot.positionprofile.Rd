\name{plot.positionprofile}
\alias{plot.positionprofile}
\alias{plot.positionprofile-method}
\alias{plot.positionprofile,cobindr-method}

\title{
function to plot a profile over the total number of predicted
transcription factor binding sites for each PWM.
}
\description{
plot.positionprofile provides position-wise profile plot over total
number of predicted TFBS for each PWM over all input
sequences. Windowing is used to provide a smoother appreance, the window
size can be adjusted with the window parameter. 
}
\usage{
\S4method{plot.positionprofile}{cobindr}(x, wind.len = 50)
}

\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences and the hits.
}
  \item{wind.len}{
    integer, defining the length of the window for counting the hits.
}
  
}

\author{
Robert Lehmann <r.lehmann@biologie.hu-berlin.de>
}
\seealso{
\code{\link{plot.positions}}
}

\keyword{ dplot }
\keyword{ hplot }
\keyword{ smooth }
