\name{pathCounts}
\alias{pathCounts}
\alias{pathCounts-method}
\alias{pathCounts,list-method}
\alias{pathCounts,procBam-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute exon path counts
}
\description{
  Compute counts for exon paths visited by aligned reads
}
\usage{
pathCounts(reads, DB, mc.cores = 1, verbose=FALSE) 
}
\arguments{
  \item{reads}{
    Object of class \code{procBam} containing aligned reads, as returned by \code{procBam}.
  }
  \item{DB}{
    Object of class \code{annotatedGenome} containing either a known or
    de novo annotated genome.
  }
  \item{mc.cores}{
    Number of processors to be used for parallel computing. Requires
    having package \code{multicore} installed and loaded.
  }
  \item{verbose}{
    Set to \code{TRUE} to print progress information.
  }
}
\section{Methods}{
  \describe{
    \item{\code{signature(reads='list')}}{
      Computes counts for exon paths from a list of procBam objects (usually reads
      processed and split by chromosome).
    }
    \item{\code{signature(reads='procBam')}}{
      Compute counts for exon paths from a procBam object of processed
      reads.
    }
  }
}
\value{
  Named integer vector with counts of exon paths. Names are character
  strings built as ".exon1.exon2-exon3.exon4.", with dashes making the
  split between exons visited by left and right-end reads
  correspondingly. 
}
\author{
  Camille Stephan-Otto Attolini
}
\seealso{
  \code{\link{procGenome}} to create an annotated genome object,
  \code{\link{createDenovoGenome}} to create a de novo annotated genome.
  See help(getNreads) to get number of fragments mapping to each island. 
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
