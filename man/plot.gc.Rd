\name{plot.gc}
\alias{plot.gc}
\alias{plot.gc-method}
\alias{plot.gc,cobindr-method}

\title{
function to visualize GC content or CpG content of input sequences
}
\description{
plot.gc calculates the GC (or CpG) content based on a window size for
each sequence and plots the content for all sequences as a heatmap over
position and sequence. 
}
\usage{
\S4method{plot.gc}{cobindr}(x, seq.ids, cpg = F, wind.size = 50,
sig.test = F, hm.margin = c(4, 10), frac = 10, n.cpu = NA)
}

\arguments{
  \item{x}{
    an object of the class "cobindr", which will hold all necessary
    information about the sequences.
}
  \item{seq.ids}{
    list of sequence identifiers, for which the GC (or CpG) content will
    be plotted. 
}
  \item{cpg}{
    logical flag, if cpg=TRUE the CpG content rather than the GC content
    will be calculated and plotted.
}
  \item{wind.size}{
    integer describing the window size for GC content calculation
}
  \item{sig.test}{
    logical flag, if sig.test=TRUE wilcoxon.test is performed per
    individual window against all windows in other sequence at the same
    position. The significance test might be slow for large number of sequences
}
  \item{hm.margin}{
    optional argument providing the margin widths for the heatmap (if sig.test=FALSE)
}
  \item{frac}{
    determines the overlap between consecutive windows as fraction wind.size/frac
  }
  \item{n.cpu}{
    number of CPUs to be used for parallelization. Default value is 'NA'
    in which case the number of available CPUs is checked and than used.
  }
}

\author{
Robert Lehmann <r.lehmann@biologie.hu-berlin.de>
}

\seealso{
\code{\link{testCpG}}
}

\examples{
library(Biostrings)

n <- 50 # number of input sequences
l <- 100 # length of sequences
bases <- c("A","C","G","T") # alphabet
# generate random input sequences with two groups with differing GC content
seqs <- sapply(1:(3*n/4), function(x) paste(sample(bases, l, replace=TRUE, 
		prob=c(.3,.22,.2,.28)), collapse=""))
seqs <- append(seqs, sapply(1:(n/4), function(x) paste(sample(bases, l, 
		replace=TRUE, prob=c(.25,.25,.25,.25)), collapse="")))
#save sample sequences in fasta file
tmp.file <- tempfile(pattern = "cobindr_sample_seq", tmpdir = tempdir(),
fileext = ".fasta")
writeXStringSet(DNAStringSet(seqs), tmp.file)

cfg <- new('configuration')
slot(cfg, 'sequence_type') <- 'fasta'
slot(cfg, 'sequence_source') <- tmp.file
# avoid complaint of validation mechanism 
slot(cfg, 'pfm_path') <- system.file('extdata/pfms',package='cobindR')
slot(cfg, 'pairs') <- ''

runObj <- new('cobindr', cfg, 'test')

plot.gc(runObj, cpg = TRUE)

unlink(tmp.file)
}

\keyword{ dplot }
\keyword{ hplot }
