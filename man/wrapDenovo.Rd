\name{wrapDenovo}
\alias{wrapDenovo}
\title{
  Run all necessary steps to get expression estimates from multiple bam files with the casper pipeline.
}
\description{
  Function to analyze bam files to generate an ExpressionSet with
  expression estimates for all samples,
  read start and fragment length distributions, path counts and
  optinally processed reads.
}
\usage{
wrapDenovo(bamFile, output_wrapKnown, knownGenomeDB, targetGenomeDB, readLength, 
  rpkm=TRUE, keep.multihits=TRUE, searchMethod="submodels", 
  exactMarginal=TRUE,  integrateMethod = "plugin", maxExons=40, 
  islandid, chroms=NULL, keep.pbam=FALSE,  keepPbamInMemory=FALSE,
  niter=10^3, priorq=3, priorqGeneExpr=2,
  mc.cores.int=1, mc.cores=1, verbose=TRUE, seed=1)
}
\arguments{
  \item{bamFile}{
    Names of bam files with the sample to analyze. These must sorted and
    indexed, and the index must be in the same directory.}
  
  \item{output_wrapKnown}{Optional argument containing the output of an
    earlier call to \code{wrapKnown}. If provided, path counts,
    read start and insert size distributions are loaded from this output
    rather than being re-computed. Better leave this argument missing
    unless you know what you're doing.}
  
  \item{knownGenomeDB}{ \code{annotatedGenome} object with known isoforms,
  e.g. from UCSC or GENCODE annotations. Used to set the prior
  probability that any given isoform is expressed. See help(calcDenovo)
  for details.}

  \item{targetGenomeDB}{\code{annotatedGenome} object with isoforms we wish to
  quantify. By default these are the same as in \code{knownGenomeDB},
  but more typically \code{targetGenomeDB} is imported from a .gtf file
  produced by some isoform prediction software.}

  \item{readLength}{Read length in bp, e.g. in a paired-end experiment where
    75bp are sequenced on each end one would set \code{readLength=75}.}

  \item{rpkm}{Set to \code{TRUE} to return reads per kilobase per million (RPKM),
    \code{FALSE} for relative expression levels.
    Important, relative expression adds up to 1 within gene island, NOT
    within gene. To get relative expressions within gene run
    \code{relexprByGene} afterwards. See help(wrapKnown).}

  \item{keep.multihits}{Set to \code{FALSE} to discard reads aligned to
  multiple positions.}

  \item{searchMethod}{Method used to perform the model search. 
  \code{"allmodels"} enumerates all possible models (warning:
  this is not feasible for genes with >5 exons). \code{"rwmcmc"} uses a
  random-walk MCMC scheme to focus on models with high posterior
  probability.
  \code{"submodels"} considers that some isoforms in \code{targetGenomeDB}
  may not be expressed, but does not search for new variants.
  \code{"auto"} uses \code{"allmodels"} for genes with up to 5
  exons and \code{"rwmcmc"} for longer genes. See help("calcDenovo").}

  \item{exactMarginal}{Set to \code{FALSE} to estimate posterior model
  probabilities as the proportion of MCMC visits. Set to \code{TRUE} to
  use the integrated likelihoods (default). See details.}

  \item{integrateMethod}{Method to compute integrated likelihoods. The default
  (\code{'plugin'}) evaluates likelihood*prior at the posterior mode and
  is the faster option. Set \code{'Laplace'} for Laplace approximations
  and \code{'IS'} for Importance Sampling. The latter increases
  computation cost very substantially.}

\item{maxExons}{Prior probabilities of isoform expression are
  estimated for genes with 1 up to \code{maxExons} exons separately,
  for genes with more than \code{maxExons} exons a combined
  estimate is used. See \code{help("modelPrior")}}

\item{islandid}{Names of the gene island to be analyzed. If missing all
  gene islands are analyzed}

\item{chroms}{Names of the chromosomes to be analyzed. If missing all
  chromosomes are analyzed.}

\item{keep.pbam}{Set to \code{TRUE} to save processed bam object, as
  returned by \code{procBam}. This object can require substantial memory
  during execution and disk storage upon saving and is not needed
  for a default analysis.}

\item{keepPbamInMemory}{Set to \code{TRUE} to keep processed bam objects
in memory to speed up some computations.}

\item{niter}{Number of MCMC iterations in the model search algorithm.}

\item{priorq}{Parameter of the Dirichlet prior for the proportion of
  reads coming from each variant. We recommend \code{priorq=3} as this
  defines a non-local prior that penalizes falsely predicted isoforms.}

\item{priorqGeneExpr}{Parameter of the Dirichlet prior distribution
  on overall gene expression. Defaults to 2 to ensure non-zero estimates.}

\item{mc.cores}{Number of cores to use in expression estimation.}

\item{mc.cores.int}{Number of cores to use when loading bam files.
  Be careful as this is a memory intensive step.}

\item{verbose}{Set to \code{TRUE} to display progress information.}
  
\item{seed}{Set seed of random number generator.}
}

\details{
  The function executes the functions \code{procBam}, \code{getDistrs},
  \code{pathCounts} \code{calcDenovo} and \code{denovoExpr}
  and formats the output nicely.
  Running \code{wrapDenovo} is much more efficient in cpu
  speed and memory usage than running these functions separately.

  When \code{rpkm} is false the function returns the estimated
  proportion of reads arising from each isoform within a gene island.
  See the details in \code{help("wrapKnown")} for more information on this.
}
\value{
  \item{denovoGenomeDB}{\code{annotatedGenome} that contains the
  isoforms in \code{targetGenomeDB} plus any new isoforms predicted by
  casper}.

  \item{exp}{Object of class \code{ExpressionSet} containing Bayesian
  model averaging expression estimates. See the fData for the posterior
  probability that each isoform is expressed.}

\item{distr}{Object of class \code{readDistrs}}

\item{pbam}{List of objects of class \code{procBam} with one element
  per chromosome}
}
\author{
  Miranda Stobbe, David Rossell
}
\references{
Rossell D, Stephan-Otto Attolini C, Kroiss M, Stocker A. Quantifying Alternative Splicing from Paired-End
RNA-sequencing data. Annals of Applied Statistics, 8(1):309-330.
}
\seealso{
\code{calcDenovo}, \code{wrapKnown}, \code{relexprByGene}
}
\examples{
## not run
## Known isoforms
##  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
##  hg19DB <- procGenome(TxDb.Hsapiens.UCSC.hg19.knownGene), genome='hg19')

## gtf with known & de novo predictions
##  mygtf <- import('hg19_denovo.gtf')
##  hg19denovoDB <- procGenome(mygtf, genome='hg19')

## bamFile="/path_to_bam/sorted.bam"
##  ans <- wrapDenovo(bamFile=bamFile, targetGenomeDB=hg19denovoDB, 
##    knownGenomeDB=hg19DB, readLength=101)

## Estimated expression via BMA
##  head(exprs(ans[['exp']]))

## Posterior probability that each isoform is expressed
##  head(fData(ans[['exp']]))
}
