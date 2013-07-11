\name{denovoGenomeExpr-class}
\Rdversion{1.1}
\docType{class}
\alias{denovoGenomeExpr-class}
\alias{show,denovoGenomeExpr-method}
\alias{[,denovoGenomeExpr,ANY,ANY,ANY-method}
\alias{[[,denovoGenomeExpr-method}
\alias{as.list,denovoGenomeExpr-method}

\title{Class "denovoGenomeExpr"}
\description{
  \code{denovoGeneExpr} stores inferred expression for de novo splicing
  variants for a single gene.
  \code{denovoGenomeExpr} stores the information for several genes
  (typically, the whole genome).
}
\section{Objects from the Class}{
  Objects are returned by \code{calcDenovo}.
}
\section{Slots}{
  \describe{
    \item{islands}{A list of \code{denovoGeneExpr} objects, with each
      element containing results for an individual gene.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "denovoGenomeExpr")}: Displays
      general information about the object. }
    \item{as.list}{Coerces the object to a list}
    \item{"["}{Selects a subset of genes}
    \item{"[["}{Selects a single gene}
	 }
}
\author{
Camille Stephan-Otto Attolini
}
\seealso{
  \code{\link{procGenome}} and \code{\link{createDenovoGenome}} to
  create \code{denovoGenomeExpr} objects.
}
\examples{
showClass("denovoGeneExpr")
showClass("denovoGenomeExpr")
}
\keyword{classes}
