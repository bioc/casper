\name{simMAE}
\alias{simMAE}
\docType{data}
\title{
  Simulate Mean Absolute Error (MAE) in estimating isoform expression
  under various experimental settings.
}
\description{
  Simulate several future RNA-seq data under various experimental
  settings (sequencing depth, read length, insert sizes), estimate
  isoform expression and assess the MAE incurred in the estimation
  process.
  The function is a wrapper combining functions simReads and calcExp.
}
\usage{
simMAE(nsim, islandid, nreads, readLength, fragLength, burnin=1000, pc, distr, readLength.pilot=readLength, eset.pilot, usePilot=FALSE, retTxsError=FALSE, genomeDB, mc.cores=1, mc.cores.int=1, verbose=FALSE, writeBam=FALSE, bamFile=NULL)
}
\arguments{
  \item{nsim}{Number of RNA-seq datasets to generate (often as little as
    \code{nsim=10} suffice)}
  \item{islandid}{When specified this argument indicates to run the
  simulations only for gene islands with identifiers in
  \code{islandid}. When not specified genome-wide simulations are
  performed.}
\item{nreads}{Vector indicating the target number of read pairs for each
  experimental setting. The actual number of reads differs from
  \code{nreads} to account for non-mappability and random read yield
  (see details)}
\item{readLength}{Vector indicating the read length in each experimental
  setting}
\item{fragLength}{Vector indicating the mean insert size in each
  experimental setting}
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
\item{writeBam}{Set to \code{TRUE} to write simulated reads to a .bam file}
\item{bamFile}{Name of the .bam file}
}
\details{
  \code{simMAE} simulates \code{nsim} datasets under each experimental
  setting defined by \code{nreads}, \code{readLength},
  \code{fragLength}. For each dataset the following steps are performed:

  1. The number of reads is nreads * readYield * pmapped, where
  readYield= runif(1,0.8,1.2) accounts for deviations in read yield and
  pmapped= runif(1,0.6,0.9)*pmappable is the proportion of mapped reads
  (60\%-90\% of the mappable reads according to the piecewise-linear
  power law of Li et al (2014))

  2. True expression levels \code{pi} are generated from their posterior
  distribution given the pilot data.

  3. Conditional on \code{pi}, RNA-seq data are generated and expression
  estimates \code{pihat} are obtained using \code{calcExp}

  4. The mean absolute estimation error \code{sum(abs(pihat-pi))} across
  all isoforms is computed


Ideally \code{simMAE} should use pilot data from a relevant related experiment
to simulate what future data may look like for the current experiment of
interest.
The recommened way to do this is to download a .bam file from such a
related experiment and processing it in casper with function
\code{wrapKnown}, as then both gene and isoform expression can be
estimated accurately.
The object output by \code{wrapKnown} is a list with elements
named \code{'pc'}, \code{'distr'} which can be given as input to
\code{simMAE}.

As an alternative to specifying \code{pc},
\code{simMAE} allows setting \code{eset.pilot} as
pilot data. Gene and isoform expression are then simulated as follows:

1. The number of reads per gene is generated from a Multinomial
distribution with success probabilities proportional to
\code{2^exprs{eset.pilot}}.

2. Relative isoform expression within each gene are generated
from a symmetric Dirichlet distribution with parameter 1/Ig, where
Ig is the number of isoforms in gene g.

We emphasize that relative isoform expressions are not trained from the
pilot data, and that while the distribution of gene expression levels
resembles that in \code{eset.pilot}, no attempt is made to match gene
identifiers and hence the results for individual genes should not be
trusted
(hence this option is only available when \code{retTxsError==FALSE}.
}
\value{
If \code{retTxsError==TRUE}, \code{simMAE} returns posterior expected MAE
  for each individual isoform.
Else the output is a \code{data.frame} with overall MAE across all isoforms
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
## maybe str(simMAE) ; plot(simMAE) ...
}
\keyword{ datagen }

