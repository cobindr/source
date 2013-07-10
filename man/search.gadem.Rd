\name{search.gadem}
\alias{search.gadem}
\alias{search.gadem-method}
\alias{search.gadem,cobindr-method}

\title{
function performs TFBS prediction denovo or based on transfac / jaspar
matrices pwms using rGADEM.
}
\description{
function performs TFBS prediction denovo or based on transfac / jaspar
matrices pwms using rGADEM. If append=T, predicted hits are appended to the hits in the input object.
}

\usage{
\S4method{search.gadem}{cobindr}(x, deNovo = FALSE, append = F, background_scan = FALSE)
}

\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences and the hits.
}
  \item{deNovo}{
    logical flag, if deNOVO=TRUE a denovo search is startet. Otherwise
  the given PFMs are used as seed.
}
  \item{append}{
    logical flag, if append=TRUE the binding sites will be appended to
  already existing results 
}
\item{background_scan}{
    logical flag, if background_scan=TRUE the function will search for
  binding sites in the set of background sequences
  }
}

\value{
  \item{x}{an object of the class "cobindr" including the predicted 
  transcription factor binding sites}
}

\references{
uses package "rGADEM" (http://www.bioconductor.org/packages/release/bioc/html/rGADEM.html)
}

\author{
Robert Lehmann <r.lehmann@biologie.hu-berlin.de>
}

\examples{
############################################################
# use simulated sequences
library(Biostrings)

n <- 600 # number of input sequences
l <- 150 # length of sequences
n.hits <- 600 # number of 'true' binding sites
bases <- c("A","C","G","T") # alphabet
# generate random input sequences with two groups with differing GC content
seqs <- sapply(1:(3*n/4), function(x) paste(sample(bases, l, replace=TRUE, 
		prob=c(.3,.22,.2,.28)), collapse=""))
seqs <- append(seqs, sapply(1:(n/4), function(x) paste(sample(bases, l, 
		replace=TRUE, prob=c(.25,.25,.25,.25)), collapse="")))
path <- system.file('extdata/pfms/myod.tfpfm',package='cobindR')
motif <- read.transfac.pfm(path)[[1]] # get PFM of binding site 
# add binding sites with distance specificity
for(position in c(70, 90)) {
	hits <- apply(apply(motif, 2, function(x) sample(x=bases, size=n.hits, 
			prob=x, replace=TRUE)), 1, paste, collapse='')
	pos.hits <- round(rnorm(n.hits, mean=position, sd=8))
	names(pos.hits) <- sample(1:n, n.hits)
	for(i in 1:n.hits) substr(seqs[as.integer(names(pos.hits)[i])], start=pos.hits[i], 
							stop=pos.hits[i]+ncol(motif)) <- hits[i] 
}
#save sample sequences in fasta file
tmp.file <- tempfile(pattern = "cobindr_sample_seq", tmpdir = tempdir(), fileext = ".fasta")
writeXStringSet(DNAStringSet(seqs), tmp.file)
#run cobindr
cfg <- cobindRConfiguration()
sequence_type(cfg) <- 'fasta'
sequence_source(cfg) <- tmp.file
sequence_origin(cfg) <- 'artificial sequences'
pfm_path(cfg) <- system.file('extdata/pfms',package='cobindR')
pairs(cfg) <- 'V$MYOD_01 V$MYOD_01' 
runObj <-cobindr(cfg, name='cobindr test using sampled sequences')

# perform tfbs prediction using rGADEM - commented out due to long time required
# runObj.bs <- search.gadem(runObj)
# show results
# plot.positions(runObj.bs)

#clean up
unlink(tmp.file)
}


\seealso{
\code{\link{rtfbs}}, \code{\link{search.pwm}}
}

\keyword{ manip }
\keyword{ iteration }
\keyword{ methods }
