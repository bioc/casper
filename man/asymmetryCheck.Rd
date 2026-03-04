\name{asymmetryCheck}
\alias{asymmetryCheck}
\alias{asymmetryCheck,data.frame-method}
\alias{asymmetryCheck,ExpressionSet-method}
\alias{asymmetryCheck,matrix-method}
\title{
  Plot asymmetry coefficients for the observed data and compare to those
  expected under Normality.
}
\description{
  Produces a boxplot for the asymmetry coefficients for each row in the
  input matrix. Normal observations are simulated using the observed sample
  means and variances, and their asymmetry coefficients are added to the plot.
}
\usage{
asymmetryCheck(x, ...)
}
\arguments{
  \item{x}{\code{ExpressionSet}, \code{matrix} or \code{data.frame}
    with genes/isoforms in rows}
  \item{...}{Other arguments to be passed on to \code{plot}}
}
\value{
  Boxplot with asymmetry coefficients for observed and simulated Normal data
}
\author{
  David Rossell
}
\examples{
mu <- rnorm(100)
x <- matrix(rnorm(100*5,mu),ncol=5)

asymmetryCheck(x)
}
\keyword{ hplot }
\keyword{distribution}
