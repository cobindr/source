\docType{package}
\name{cobindR-package}
\alias{cobindr}
\alias{cobindr-package}
\alias{cobindR}
\alias{cobindR-package}
\title{An R package for analyzing co-occurring transcription factor binding sites}
\description{
Many transcription factors (TFs) regulate gene expression by binding to
specific DNA motifs near genes. Often the regulation of gene expression
is not only controlled by one TF, but by many TFs together, that can
either interact in a cooperative manner or interfere with each other. In
recent years high thoughput methods, like ChIP-Seq, have become
available to produce large amounts of data, that contain potential
regulatory regions. In silico analysis of trancription factor binding
sites can help to interpret these enormous datasets in a convenient and
fast way or narrow down the results to the most significant regions for
further experimental studies.

cobindR provides a complete set of methods to analyse and detect pairs
of TFs, including support of diverse input formats and different
background models for statistical testing. Several visualization tools
are implemented to ease the interpretation of the results. 
}

\author{
  Yue-Hien Lee,  Robert Lehmann,  Stefan Kroeger,  Manuela Benary 
}

\seealso{
  The core class in this package: \code{\link{cobindr-class}}. 
  The core function in this package: \code{\link{find.pairs}}. 
}
