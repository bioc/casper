\name{getDistrs}
\alias{getDistrs}
\title{
  Compute fragment start and fragment length distributions
}
\description{
  Compute fragment start distributions by using reads
  aligned to genes with only one annotated variant.
  Estimate fragment length
  distribution using fragments aligned to long exons (>1000nt).
  Fragment length is defined as the distance between the start of the
  left-end read and the end of the right-end read.
}
\usage{
getDistrs(DB, bam, pbam, islandid, verbose, nreads=4*10^6, readLength)
}
\arguments{
  \item{DB}{
    Annotated genome. Object of class \code{knownGenome} as returned by \code{procGenome}.
  }
  \item{bam}{
    Aligned reads, as returned by \code{scanBam}. It must be a list with
    elements 'qname', 'rname', 'pos' and 'mpos'. Ignored when argument \code{pbam} is
    specified.
  }
  \item{pbam}{
    Processed BAM object of class \code{procBam}, as returned by
    function \code{procBam}. Arguments \code{bam} and \code{readLength}
    are ignored when \code{pbam} is specified. 
  }
  \item{islandid}{
    Island IDs of islands to be used in the read start distribution
    calculations (defaults to genes with only one annotated variant)
  }
  \item{verbose}{
    Set to \code{TRUE} to print progress information.
  }
  \item{nreads}{
    To speed up computations, only the first \code{nreads} are used to
    obtain the estimates. The default value of 4 milions usually gives
    highly precise estimates.
  }
  \item{readLength}{
    Read length in bp, e.g. in a paired-end experiment where
  75bp are sequenced on each end one would set \code{readLength=75}.
  }
}
\value{
  An object of class readDistrs with slots:
  \item{lenDis}{Table with number of fragments with a given length}
  \item{stDis}{Cumulative distribution function (object of type closure) for relative start position}
}
\author{
Camille Stephan-Otto Attolini, David Rossell
}
\examples{
data(K562.r1l1)
data(hg19DB)
bam0 <- rmShortInserts(K562.r1l1, isizeMin=100)

distrs <- getDistrs(hg19DB,bam=bam0,readLength=75)

#Fragment length distribution
plot(distrs,'fragLength')

#Fragment start distribution (relative to transcript length)
plot(distrs,'readSt')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stats }

