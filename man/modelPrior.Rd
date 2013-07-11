\name{modelPrior}
\alias{modelPrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set prior distribution on expressed splicing variants.
}
\description{
  Set prior on expressed splicing variants using the genome annotation
  contained in a \code{knownGenome} object.

  The prior probability of variants V1,...,Vn being expressed depends on
  n, on the number of exons in each variant V1,...,Vn and 
  the number of exons in the gene.
  See the details section.
}
\usage{
modelPrior(genomeDB, maxExons=40, smooth=TRUE, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genomeDB}{Object of class \code{knownGenome} }
  \item{maxExons}{The prior distribution is estimated for genes with 1
  up to \code{maxExons} exons. As there are fewer genes with many exons,
  the prior parameters are estimated poorly. To avoid this common
  estimate is used for all genes with more than \code{maxExons} exons}
  \item{smooth}{If set to \code{TRUE} the estimated prior distribution
  parameters for the number of exons in a gene are smoothed using
  Generalized Additive Models. This step typically improves the
  precision of the estimates, and is only applied to genes with 10 or
  more exons.}
  \item{verbose}{Set to \code{TRUE} to print progress information.}
}
\details{
  The goal is to set a prior that takes into account the number of
  annotated variants for genes with E exons, as well as the number of exons in
  each variant.
  
  Suppose we have a gene with E exons.
  Let V_1,...,V_n be n variants of interest and let |V_1|,...,|V_n| be the
  corresponding number of exons in each variant.
  The prior probability of variants V_1,...,V_n being expressed is modeled
  as

  P(V_1,...,V_n|E)= P(n|E) P(|V_1| |E) ... P(|V_n| |E)

  where P(n|E)= NegBinom(n; k_E, r_E) I(0 < n < 2^E) and
  P(|V_i| |E)= BetaBinomial(|V_i|-1; E-1, alpha_E, beta_E).
  
  The parameters k_E, r_E, alpha_E, beta_E depend on E (the number of exons
  in the gene) and are estimated from the available annotation via
  maximum likelihood.
  Parameters are estimated jointly for all genes with E>=
  \code{maxExons} in order to improve the precision.
  
  For \code{smooth==TRUE}, alpha_E and beta_E are modeled as a smooth
  function of E by calling \code{gam} and setting the smoothing
  parameter via cross-validation. Estimates for genes with E>=10 are
  substituted by their smooth versions, which typically helps improve
  stability in the estimates.
}
\value{
  List with 2 components.
  
  \item{nvarPrior}{List with prior distribution on the number of
  expressed variants for genes with 1,2,3... exons. Each element contains
  the truncated Negative Binomial parameters, observed and predicted
  frequencies (counting the number of genes with a given number of
  variants).}

  \item{nexonPrior}{List with prior distribution on the number of exons
  in a variant for genes with 1,2,3... exons. Each element contains the
  Beta-Binomial parameters, observed and predicted frequencies (counting
  the number of variants with a given number of exons)}
}
\author{
David Rossell, Camille Stephan-Otto Attolini
}
\examples{
data(hg19DB)
mprior <- modelPrior(hg19DB, maxExons=10)

##Prior on number of expressed variants
##Genes with 2 exons
##mprior$nvarPrior[['2']]
##Genes with 3 exons
##mprior$nvarPrior[['3']]

##Prior on the number of exons in an expressed variant
##Genes with 2 exons
##mprior$nexonPrior[['2']]
##Genes with 3 exons
##mprior$nexonPrior[['3']]
}
\keyword{ stats }
