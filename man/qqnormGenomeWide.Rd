\name{qqnormGenomeWide}
\alias{qqnormGenomeWide}
\alias{qqnormGenomeWide,data.frame-method}
\alias{qqnormGenomeWide,ExpressionSet-method}
\alias{qqnormGenomeWide,matrix-method}
\alias{qqgammaGenomeWide}
\alias{qqgammaGenomeWide,data.frame-method}
\alias{qqgammaGenomeWide,ExpressionSet-method}
\alias{qqgammaGenomeWide,matrix-method}
\title{
  Genome-wide qq-normal and qq-gamma plots
}
\description{
  qqnormGenomeWide overlays quantile-quantile normal plots (qqnorm) for a series of genes
  (rows in the input matrix), to provide an overall assessment of
  Normality. Similarly, qqgammaGenomeWide overlays quantile-quantile gamma
  plots.

  Note that the theoretical quantiles for z-scores under
  a Normal are the same for all genes, but the gamma theoretical quantiles
  depend on the Gamma parameter estimates for each gene and hence the
  theoretical quantiles are different for each gene (resulting in
  different x-values in each qq-plot)
}
\usage{
qqnormGenomeWide(x, ngenes=min(1000, nrow(x)), ...)

qqgammaGenomeWide(x, ngenes=min(1000, nrow(x)), ...)
}
\arguments{
  \item{x}{\code{ExpressionSet}, \code{matrix} or \code{data.frame}
    with genes/isoforms in rows}
  \item{ngenes}{A qqnorm plot is produced for the first \code{ngenes}
    rows in \code{x}}
  \item{...}{Other arguments to be passed on to \code{plot}}
}
\value{
  Produces a figure overlaying qq-normal or qq-gamma plots for \code{ngenes}
  comparing observed vs. theoretical quantiles
}
\author{
  David Rossell
}
\examples{
mu <- rnorm(100)
x <- matrix(rnorm(100*5,mu),ncol=5)

qqnormGenomeWide(x)
qqgammaGenomeWide(exp(x))
}
\keyword{ hplot }
\keyword{distribution}
