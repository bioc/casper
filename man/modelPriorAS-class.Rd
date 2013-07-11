\name{modelPriorAS-class}
\Rdversion{1.1}
\docType{class}
\alias{modelPriorAS-class}
\alias{show,modelPriorAS-method}
\alias{[,modelPriorAS,ANY,ANY,ANY-method}
\alias{coef,modelPriorAS,ANY,ANY,ANY-method}

\title{Class "modelPriorAS"}
\description{
  \code{modelPriorAS} stores parameters for the prior distribution on
  all possible alternative splicing configuration (i.e. prior on model
  space). This information is used for de novo reconstruction of
  splicing variants.
}
\section{Objects from the Class}{
  Objects are created by function \code{modelPrior}.
}
\section{Slots}{
  \describe{
    \item{nvarPrior}{Prior on the number of variants per gene. A list with components \code{"nbpar"} containing
      the parameters of the Negative Binomial distribution, \code{"obs"}
      containing the observed counts and \code{"pred"} the Negative
      Binomial predicted counts.}
    \item{nexonPrior}{Prior on the number of exons in an expressed
      variant. A list with components \code{"bbpar"} containing
      Beta-Binomial parameters, \code{"obs"} containing the observed
      counts and \code{"pred"} the Beta-Binomial predicted counts.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "modelPriorAS")}: Displays
      general information about the object. }
    \item{"["}{Selects prior parameters for genes with the specified
      number of exons}
    \item{coef}{Selects a single gene}
	 }
}
\author{
David Rossell
}
\seealso{
  \code{\link{procGenome}} and \code{\link{createDenovoGenome}} to
  create \code{modelPriorAS} objects.
}
\examples{
showClass("modelPriorAS")
}
\keyword{classes}
