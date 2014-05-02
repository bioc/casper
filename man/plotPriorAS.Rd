\name{plotPriorAS}
\alias{plotPriorAS}
\alias{plotPriorAS-methods}
\alias{plotPriorAS,modelPriorAS-method}
\docType{methods}
\title{
Plot prior distribution on set of expressed variants (i.e. the model space).
}
\description{
Plots the prior distribution on the number of expressed variants and the
number of exons per variant in genes with \code{exons} exons
(as returned by function \code{modelPrior}).
The prior distribution is compared to the observed frequencies to check
that the assumed distributional forms are reasonable.
}
\usage{
plotPriorAS(object, type="nbVariants", exons=1:9, xlab,
ylab="Probability", col=c("red","blue"))
}
\arguments{
 \item{object}{\code{modelPriorAS} object with prior distribution on
   model space.}
 \item{type}{Set to \code{"nbVariants"} to plot the prior on the number
   of variants per gene. Set to \code{"nbExons"} to plot the prior on
   the number of exons.}
 \item{exons}{Vector with integers. The plot is only produced with
   number of exons indicated in \code{exons}.}
 \item{xlab}{x-axis label, passed on to \code{plot}}
 \item{ylab}{y-axis label, passed on to \code{plot}}
 \item{col}{Colors for bars showing prior probabilities and frequencies
   in the known genome}
}
\section{Methods}{
\describe{
\item{\code{signature(object = "modelPriorAS")}}{
  \code{object} contains the prior distribution on the model space, as
  returned by function \code{modelPrior}
  }
}}
\examples{
#See modelPrior examples
}
\keyword{hplots}
