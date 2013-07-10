\name{testCpG}
\alias{testCpG}
\alias{testCpG-method}
\alias{testCpG,cobindr-method}

\title{
function to cluster sequences based on their CpG and GC content
}
\description{
diagnostical function - GC content and CpG content are clustered using 2D gaussian
 models (Mclust). FALSE is returned if > max.clust (default=1) subgroups are found
 using the bayesian information criterion (BIC). If do.plot=TRUE, the results are visualized.
}
\usage{
\S4method{testCpG}{cobindr}(x, max.clust = 4, do.plot = F, n.cpu = NA)
}

\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences and the hits.
}
  \item{max.clust}{
    integer describing the maximal number of clusters which are used for
    separating the data. 
}
  \item{do.plot}{
    logical flag, if do.plot=TRUE a scatterplot for the GC and CpG
    content for each sequence is produced and the clusters are color coded.
  }
  \item{n.cpu}{
    number of CPUs to be used for parallelization. Default value is 'NA'
    in which case the number of available CPUs is checked and than used.
  }
}

\value{
 \item{result}{logical flag, FALSE is returned if more than one subgroups are found
 using the bayesian information criterion (BIC)}
 \item{gc}{matrix with rows corresponding to sequences and columns
 corresponding to GC and CpG content}
}
\references{
the method uses clustering functions from the package "mclust" (http://www.stat.washington.edu/mclust/)
}
\author{
Robert Lehmann <r.lehmann@biologie.hu-berlin.de>
}

\seealso{
\code{\link{plot.gc}}
}

\examples{
cfg <- cobindRConfiguration()
sequence_type(cfg) <- 'fasta'
sequence_source(cfg) <- system.file('extdata/example.fasta', package='cobindR')
# avoid complaint of validation mechanism 
pfm_path(cfg) <- system.file('extdata/pfms',package='cobindR')
pairs(cfg) <- '' 
runObj <- cobindr( cfg)
testCpG(runObj, max.clust = 2, do.plot = TRUE) 
}

\keyword{dplot}
\keyword{manip}
\keyword{clust}
