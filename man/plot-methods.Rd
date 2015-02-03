\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot,readDistrs,ANY-method}
\alias{lines}
\alias{lines,readDistrs-method}
\title{
Plot estimated read start and fragment length distributions.
}
\description{
  Plots the estimated fragment length (insert size) distribution
  and the relative read start distribution
  (0 indicating transcription start, 1 transcription end).
  The former checks that the insert size distribution
  matches that described in the experimental protocol.
  The latter checks the extent to which reads are non-uniformly
  distributed (note: casper does NOT assume reads to be uniformly
  distributed, so a lack of uniformity is not a problem per se).
}
\arguments{
 \item{x}{Object of type \code{readDistrs}, as returned by \code{getDistrs}.}
 \item{y}{Set to \code{"fragLength"} to plot the estimated insert size
   ditribution. Set to \code{"readSt"} to plot a histogram of the
   estimated read start distribution.}
 \item{...}{Further arguments to be passed on to plot.}
}
\section{Methods}{
  \describe{

  \item{\code{signature(x = "readDistrs", y = "ANY")}}{
    \code{x} is an object of type \code{readDistrs}, as returned by
    \code{getDistrs}. The plot function allows to visualize the fragment
    length and read start distributions.
  }

  \item{\code{signature(x = "readDistrs")}}{
    \code{x} is an object of type \code{readDistrs}, as returned by
    \code{getDistrs}. The plot function allows to visualize the fragment
    length and read start distributions.
  }
}}
\examples{
#See getDistrs examples
}
\keyword{hplots}
