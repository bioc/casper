\name{simulatedSamples-class}
\Rdversion{1.1}
\docType{class}
\alias{simulatedSamples-class}
\alias{show,simulatedSamples-method}
\alias{coef,simulatedSamples-method}
\alias{exprs,simulatedSamples-method}
\alias{[,simulatedSamples,ANY,ANY,ANY-method}
\alias{[,simulatedSamples,missing,index,missing-method}


\title{Class "simulatedSamples"}
\description{
  \code{simulatedSamples} stores multiple simulated isoform expression
  datasets. Each dataset contains the (simulation) true mean expression
  in each group and residual variance, as well as the estimated
  expression in each individual sample.
}
\section{Objects from the Class}{
  Objects are returned by \code{simMultSamples}.
}
\section{Slots}{
  The class extends a list directly.
  \describe{
    \item{.Data}{A list, each element containing a different simulated dataset}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "simulatedSamples")}: Displays
      general information about the object. }
    \item{coef}{\code{signature(object = "simulatedSamples")}:
      Returns a matrix with difference between group means (simulation
      truth) in all simulated datasets}
    \item{exprs}{\code{signature(object = "simulatedSamples")}: Returns
      a list of ExpressionSets containing the estimated expressions in
      each simulation.}
    \item{mergeBatches}{\code{signature(x="ExpressionSet",y="simulatedSamples")}:
    Combines \code{x} with each element in \code{exprs} in \code{y}, and
    returns a list. See \code{help(mergeBatches)} for more details.}
    \item{"["}{\code{x[i]} selects a subset of simulations, \code{x[,j]}
    a subset of the samples in each simulation}
  }
}
\author{
David Rossell
}
\seealso{
  \code{\link{mergeBatches}}
}
\examples{
showClass("simulatedSamples")
}
\keyword{classes}
