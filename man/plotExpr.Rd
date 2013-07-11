\name{plotExpr}
\alias{plotExpr}
\alias{plotExpr-methods}
\alias{plotExpr,denovoGeneExpr-method}
\docType{methods}
\title{
Plot inferred gene structure and expression.
}
\description{
Plots variants with sufficiently large posterior probability of being
expressed along with their (marginal) estimated expression.
}
\usage{
plotExpr(gene, minProbExpr = 0.5, minExpr = 0.1,
    xlab = "(kb)", ylab = "", xlim, cex = 1, yaxt = "n", col, ...) 
}
\arguments{
 \item{gene}{\code{denovoGeneExpr} object containing results for a
   single gene, as returned by
   \code{calcDenovo}.}
 \item{minProbExpr}{Variants with marginal posterior probability of
   expression below \code{minProbExpr} are not reported}
 \item{minExpr}{Variants with (marginal) estimated expression below
   \code{minExpr} are not reported. Can be useful to remove sequence
   preference artifacts.}
 \item{xlab}{x-axis label, passed on to \code{plot}}
 \item{ylab}{y-axis label, passed on to \code{plot}}
 \item{xlim}{x-axis limits, passed on to \code{plot}}
 \item{cex}{Character expansion, passed on to \code{plot}}
 \item{yaxt}{Type of y-axis, passed on to \code{plot}}
 \item{col}{Colors for each variant, defaults to rainbow colors. It is
   possible to specify a single color.}
 \item{...}{Other arguments to be passed on to \code{plot}}
}
\section{Methods}{
\describe{
\item{\code{signature(gene = "denovoGeneExpr")}}{
  \code{gene} contains the results from a de novo isoform expression
  analysis for a single gene, as returned by \code{calcDenovo}. When
  \code{calcDenovo} is run on multiple genes simultaneously, the
  desired gene can be selected using the "[[" operator as usual.
  }
}}

\details{
  The marginal posterior probability that a variant is expressed is the
  sum of the posterior probabilities of all models containing that
  variant.
  
  The marginal estimated expression is the average expression across all
  models (including those where the variant has 0 expression) weighted
  by the posterior probability of each model.
}
\examples{
#See calcDenovo examples
}
\keyword{hplots}
