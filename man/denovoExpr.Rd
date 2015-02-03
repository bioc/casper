\name{denovoExpr}
\alias{denovoExpr}
\title{
  Estimate expression for de novo splicing variants. 
}
\description{
  Obtains expression estimates from \code{denovoGenomeExpr} objects, as
  returned by \code{calcDenovo}.
  When \code{rpkm} is set to \code{TRUE}, fragments per kilobase
  per million are returned. Otherwise relative expression estimates are
  returned.

  The estimates can be obtained by Bayesian model averaging (default)
  or by selecting the model with highest posterior probability. See details.
}
\usage{
denovoExpr(x, pc, rpkm = TRUE, summarize = "modelAvg", minProbExpr = 0.5, minExpr = 0.05) 
}
\arguments{
\item{x}{ \code{denovoGenomeExpr} object returned by \code{calcExp} }
\item{pc}{Named vector of exon path counts as returned by \code{pathCounts}}
\item{rpkm}{Set to \code{FALSE} to return relative expression levels,
  i.e. the proportion of reads generated from each variant per
  gene. These proportions add up to 1 for each gene.
  Set to \code{TRUE} to return fragments per kilobase per million
  (RPKM).}
\item{summarize}{Set to \code{"modelAvg"} to obtain model averaging
  estimates, or to \code{"bestModel"} to select the model with highest
  posterior probability. We recommend the former, as even the best model
  may have low posterior probability.}
\item{minProbExpr}{Variants with (marginal posterior) probability of being
  expressed below \code{minProbExpr} are omitted from the results. This
  argument is useful to eliminate variants that are not at least
  moderately supported by the data.}
\item{minExpr}{Variants with relative expression \code{minExpr} are
  omitted from the results. This is useful to eliminate variants to
  which few reads are assigned, e.g. due to read miss-alignments or biases.}
}
\value{
  Expression set with expression estimates.
  The \code{featureData} indicates the gene island id, posterior probability
  that each variant is expressed (column \code{"probExpressed"}), the
  number of aligned reads per gene island (column \code{"explCnts"}) and the exons in each variant.
}
\references{
Rossell D, Stephan-Otto Attolini C, Kroiss M, Stocker A. Quantifying Alternative Splicing from Paired-End
RNA-sequencing data. COBRA (http://biostats.bepress.com/cobra/art97)
}
\author{
  David Rossell
}
\examples{

## NOTE: toy example with few reads & genes to illustrate code usage
##       Results with complete data are much more interesting!

data(K562.r1l1)
data(hg19DB)

#Pre-process
bam0 <- rmShortInserts(K562.r1l1, isizeMin=100)
pbam0 <- procBam(bam0)

#Estimate distributions, get path counts
distrs <- getDistrs(hg19DB,bam=bam0,readLength=75)
pc <- pathCounts(pbam0, DB=hg19DB)

#Set prior distrib on model space
mprior <- modelPrior(hg19DB, maxExons=40, smooth=FALSE)

#Fit model
denovo <- calcDenovo(distrs,targetGenomeDB=hg19DB,pc=pc,readLength=75,priorq=3,mprior=mprior,minpp=0)

head(names(denovo))
denovo[['6499']]
posprob(denovo[['6499']])

#Get estimates
eset <- denovoExpr(denovo, pc=pc, rpkm=TRUE, minProbExpr=0.5)

head(exprs(eset))
head(fData(eset))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ htest }% __ONLY ONE__ keyword per line
