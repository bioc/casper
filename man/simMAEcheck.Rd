\name{simMAEcheck}
\alias{simMAEcheck}
\docType{data}
\title{
  Model checking for One Sample Problems.
}
\description{
  Simulates RNA-seq data under the same experimental setting as in the observed data,
  and compares the observed vector of number of reads per gene with the simulations.
}
\usage{simMAEcheck(nsim, islandid, burnin=1000, pc, distr, readLength.pilot, 
eset.pilot, usePilot=FALSE, retTxsError=FALSE, genomeDB, mc.cores=1, 
mc.cores.int=1, verbose=FALSE)}
\arguments{
  \item{nsim}{Number of RNA-seq datasets to generate (often as little as
    \code{nsim=10} suffice)}
  \item{islandid}{When specified this argument indicates to run the
  simulations only for gene islands with identifiers in
  \code{islandid}. When not specified genome-wide simulations are
  performed.}
\item{burnin}{Number of MCMC burn-in samples (passed on to \code{calcExp})}
\item{pc}{Observed path counts in pilot data. When not specified, these
  are simulated from \code{eset.pilot}}
\item{distr}{Estimated read start and insert size distributions in pilot
  data}
\item{readLength.pilot}{Read length in pilot data}
\item{eset.pilot}{ExpressionSet with pilot data expression in
  log2-RPKM, used to simulate \code{pc} when not specified by the
  user. See details}
\item{usePilot}{By default \code{casper} assumes that the pilot data is from a related experiment rather than the current tissue of interest (\code{usePilot=FALSE}). Hence, the pilot data is used to simulate new RNA-seq data but not to estimate its expression. However, in some cases we may be interested in re-sequencing the pilot sample at deeper length, in which case one would want to combine the pilot data with the new data to obtain more precise estimates. This can be achieved by setting \code{usePilot=TRUE}}
\item{retTxsError}{If \code{retTxsError=TRUE}, \code{simMAE} returns posterior expected MAE
for each individual isoform. This option is not available when
\code{eset.pilot} is specified instead of \code{pc}.
Else the output is a \code{data.frame} with overall MAE across all isoforms}
\item{genomeDB}{\code{annotatedGenome} object, as returned by
  \code{procGenome}}
\item{mc.cores}{Number of cores to use in the expression estimation step, passed on to \code{calcExp}}
\item{mc.cores.int}{Number of cores to simulate RNA-seq datasets in parallel}
\item{verbose}{Set \code{verbose=TRUE} to print progress information}
}
\details{
  \code{simMAEcheck} simulates \code{nsim} datasets under the same experimental setting
  as in the observed data. For more details, please check the documentation for 
  \code{simMAE}, which is the basis of this function.
}
\value{
The output is a list with 2 entries. The first entry is a \code{data.frame} 
with overall MAE across all isoforms in the simulations (see \code{simMAE} for details). 
The second entry contains the expected number of genes for which the number of
reads in the data lies in the range of the posterior predictive simulations (under the hypothesis that they have the same
distribution) and the actual number of genes for which the condition is satisfied.
}
\seealso{wrapKnown,simReads,calcExp}
\references{
  Stephan-Otto Attolini C., Pena V., Rossell D. Bayesian designs for
  personalized alternative splicing RNA-seq studies (2014)

  Li, W. and Freudenberg, J. and Miramontes, P. Diminishing return
  for increased Mappability with longer sequencing reads: implications of
  the k-mer distributions in the human genome. BMC Bioinformatics, 15,
  2 (2014)
}
\examples{
#Run casperDesign() to see full manual with examples
}
\keyword{ datagen }

