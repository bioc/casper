\name{splitGenomeByLength}
\alias{splitGenomeByLength}
\title{
  Split an annotatedGenome object into subsets according to gene length
}
\description{
  \code{splitGenomeByLength} splits an annotatedGenome according to gene
  length (bp), which allows estimating the fragment start and length
  distribution for each subset separately.
}
\usage{
splitGenomeByLength(DB, breaks=c(0,3000,5000,Inf))
}
\arguments{
  \item{DB}{Object containing annotated genome. Must be of class
    \code{annotatedGenome}, as returned by \code{procGenome} or \code{createDenovoGenome}.}
  \item{breaks}{Breakpoints to define gene subgroups.}
}
\details{
   By default groups are <3000bp, 3000-5000bp, >5000bp, which work well for the human
    genome. Further sub-dividisions may result in unstable estimates of
    fragment start and length distributions.
}
\value{
  List where each component is of class \code{annotatedGenome}.
}
\author{
  David Rossell
}
\seealso{
  \code{procGenome} and \code{createDenovoGenome} for creating
  \code{annotatedGenome} objects.
  \code{getDistrs} for estimating fragment start and length distribution.
}
\examples{
##Not run
## genDB<-makeTranscriptDbFromUCSC(genome="hg19", tablename="refGene")
## hg19DB <- procGenome(genDB, "hg19")
## hg19split <- splitGenomeByLength(hg19DB)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~annotation }
