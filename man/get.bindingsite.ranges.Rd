\name{get.bindingsite.ranges}
\alias{get.bindingsite.ranges}

\title{
convenience function to convert predicted binding sites to GRanges object.
}
\description{
Function converts predicted binding sites into a GRanges object (package: GenomicFeatures). This allows for
easy interaction with other tools as well as output of different formats (bed, gff).
}
\usage{
get.bindingsite.ranges(x, ...)
}

\arguments{
  \item{x}{
	An object of the class "cobindr", which will hold the predicted binding site locations.
      }
  \item{...}{optional additional parameters}    
}

\value{
	A GRanges object holding the positions of all predicted transcription factor binding sites relative to 
	the input sequence.

}

\author{
Robert Lehmann <r.lehmann@biologie.hu-berlin.de>
}



\seealso{
get.pairs
write.bindingsites
write.bindingsites.table
}

\examples{
# export(get.bindingsite.ranges(runObj), "tfbs_hits.gff3")
}
\keyword{ manip }
\keyword{ methods }
