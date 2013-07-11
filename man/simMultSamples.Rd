\name{simMultSamples}
\alias{simMultSamples}
\title{
  Simulate paired end reads for multiple samples
  based on pilot data
}
\description{
  Draw true expression levels for each group and future individual
  samples within each group, and simulate future observed data (in the
  form of Casper expression estimates).
  
  These simulations can serve as the basis for sample size calculation:
  if one were to sequence \code{nsamples} new RNA-seq samples, what data
  would we expect to see? The simulation is posterior predictive,
  i.e. based on the current available data \code{x}.
}
\usage{
simMultSamples(B, nsamples, nreads, readLength, x, groups='group', distrs, genomeDB, verbose=TRUE, mc.cores=1)
}
\arguments{
  \item{B}{Number of simulations to obtain}
  \item{nsamples}{Vector indicating number of future samples per group,
    e.g. \code{nsamples=c(5,5)} to simulate 5 new samples for 2 groups.}
  \item{nreads}{number of RNA-seq reads per
  sample. \code{simMultSamples} assumes that exactly these many reads will be
  aligned to the genome.}
  \item{readLength}{Read length, i.e. in an experiment with paired reads
    at 100bp each, \code{readLength=100}.}
  \item{x}{\code{ExpressionSet} containing pilot data. \code{x[[group]]} indicates groups to be compared}
  \item{groups}{Name of column in \code{pData(x)} indicating the groups}
  \item{distrs}{Fragment start and length distributions. It can be
  either an object of class readDistrs, or a list where each element is
  of class readDistrs. 
  In the latter case, an element is chosen at random for each individual sample (so that uncertainty in these distributions is taken into account)}
  \item{genomeDB}{annotatedGenome object}
  \item{verbose}{Set to \code{TRUE} to print progress}
  \item{mc.cores}{Number of cores to use in function.
  \code{mc.cores>1} requires package \code{multicore}}
}

\value{
  List of length \code{B}, where each element is itself a list containing an individual simulation.
  
  \item{simTruth}{\code{data.frame} indicating the mean and
  standard deviation of the Normal distribution used to generate data
  from each group}
  \item{simExpr}{\code{ExpressionSet} with Casper expression estimates,
  as returned by \code{calcExp}. \code{pData(simExpr)} indicates group
  information, and \code{fData(simExpr)} the number of observed reads
  for each gene}
}
\author{
Victor Pena
}
\keyword{ datagen }
