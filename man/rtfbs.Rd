\name{rtfbs}
\alias{rtfbs}
\alias{rtfbs-method}
\alias{rtfbs,cobindr-method}

\title{
function performs TFBS prediction using the package rtfbs
}
\description{
function performs TFBS prediction using the package rtfbs
}
\usage{
\S4method{rtfbs}{cobindr}(x, append = F, background_scan = FALSE, n.cpu = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences and the hits.
}
  \item{append}{
    logical flag, if append=TRUE the binding sites will be appended to
  already existing results 
}
  \item{background_scan}{
    logical flag, if background_scan=TRUE the background sequences will
  be searched for transcription factor binding sites
}
  \item{n.cpu}{
    number of CPUs to be used for parallelization. Default value is 'NA'
    in which case the number of available CPUs is checked and than used.
  }
}

\value{
  \item{x}{an object of the class "cobindr" including the predicted 
  transcription factor binding sites}
}

\references{
uses the package "rtfbs" (http://cran.r-project.org/web/packages/rtfbs/index.html)
}

\author{
Yue-Hien Lee <>
}

\examples{

############################################################
# use simulated sequences
library(Biostrings)

n <- 400 # number of input sequences
l <- 500 # length of sequences
n.hits <- 250 # number of 'true' binding sites
bases <- c("A","C","G","T") # alphabet
# generate random input sequences with two groups with differing GC content
seqs <- sapply(1:(3*n/4), function(x) paste(sample(bases, l, replace=TRUE, 
		prob=c(.3,.22,.2,.28)), collapse=""))
seqs <- append(seqs, sapply(1:(n/4), function(x) paste(sample(bases, l, 
		replace=TRUE, prob=c(.25,.25,.25,.25)), collapse="")))
path <- system.file('extdata/pfms/myod.tfpfm',package='cobindR')
motif <- read.transfac.pfm(path)[[1]] # get PFM of binding site 
# add binding sites with distance specificity
for(position in c(110, 150)) {
	hits <- apply(apply(motif, 2, function(x) sample(x=bases, size=n.hits, 
			prob=x, replace=TRUE)), 1, paste, collapse='')
	pos.hits <- round(rnorm(n.hits, mean=position, sd=8))
	names(pos.hits) <- sample(1:n, n.hits)
	for(i in 1:n.hits) substr(seqs[as.integer(names(pos.hits)[i])], 
						start=pos.hits[i], stop=pos.hits[i]+ncol(motif)) <- hits[i] 
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
fdrThreshold(cfg) <- 0 
runObj <- cobindr(cfg, name='cobindr test using sampled sequences')
# perform tfbs prediction using rtfbs
runObj.bs <- rtfbs(runObj)
# show results
plot.positionprofile(runObj.bs)

#clean up
unlink(tmp.file)
}

\seealso{
\code{\link{search.pwm}}, \code{\link{search.gadem}}
}

\keyword{ manip }
\keyword{ iteration }
\keyword{ methods }

