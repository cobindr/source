\name{plot.tfbslogo}
\alias{plot.tfbslogo}
\alias{plot.tfbslogo-method}
\alias{plot.tfbslogo,cobindr-method}

\title{
function to plot sequence logos based on hits of tools
}

\description{
plot.tfbslogo produces a sequence logo based on all hits per position
weight matrix. If a file path is specified in pdf.name, sequences logos
will be written into the specified file.
}

\usage{
\S4method{plot.tfbslogo}{cobindr}(x, pwms)
}

\arguments{
  \item{x}{Object}
  \item{pwms}{
    vector of names of position weight matrices used for searching the
    sequences. For each pwm a new sequence logo based on the hits is produced.
}
}

\author{
Robert Lehmann <r.lehmann@biologie.hu-berlin.de
}



