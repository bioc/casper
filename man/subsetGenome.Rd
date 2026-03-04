\name{subsetGenome}
\alias{subsetGenome}
\alias{subsetGenome-methods}
\alias{subsetGenome,character,missing,annotatedGenome-method}
\alias{subsetGenome,missing,character,annotatedGenome-method}
\docType{methods}
\title{
  subsetGenome subsets an object of class annotatedGenome for a set of
  island IDs or chromosome names.
}
\description{
  ~~ Methods for function \code{subsetGenome} in package \pkg{casper} ~~
  Subset an annotatedGenome object by islands or chromosomes.
  
}
\usage{
subsetGenome(islands, chr, genomeDB)
}
\arguments{
  \item{islands}{Vector of characters with the island IDs to retrieve
    from genome.}
  \item{chr}{Vector of characters with the names of chromosomes to
    retrieve from genome.}
  \item{genomeDB}{annotatedGenome object with genome to subset.}
}
\section{Methods}{
\describe{

\item{\code{signature(islands = "character", chr = "missing", genomeDB = "annotatedGenome")}}{
  %%  ~~describe this method here~~
  Subset annotatedGenome object by a set of island IDs.
}

\item{\code{signature(islands = "missing", chr = "character", genomeDB = "annotatedGenome")}}{
  %%  ~~describe this method here~~
  Subset annotatedGenome object by chromosomes.
}
}}
\keyword{methods}

