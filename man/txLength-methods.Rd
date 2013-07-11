\name{txLength}
\alias{txLength}
\alias{txLength-methods}
\alias{txLength,character,missing,annotatedGenome-method}
\alias{txLength,missing,character,annotatedGenome-method}
\alias{txLength,missing,missing,annotatedGenome-method}
\docType{methods}
\title{ ~~ Methods for Function \code{txLength} in Package \pkg{casper} ~~}
\description{
 ~~ Methods for function \code{txLength} in package \pkg{casper} ~~
 Function to retrieve transcript lengths from annotated genome (class genomeDB).
}

\usage{
  txLength(islandid, txid, genomeDB)
}
\arguments{
  \item{islandid}{
    Retrieve length for transcripts in island \code{islandid}.
  }
  \item{txid}{
    Retrieve length for \code{txid} transcripts.
  }
  \item{genomeDB}{
    Annotated genome of class \code{genomeDB}.
    }
  }
  \section{Methods}{
    \describe{
      
      \item{\code{signature(islandid = "character", txid = "missing", genomeDB = "annotatedGenome")}}{
	%%  ~~describe this method here~~
	Retrieve lengths from genomeDB for transcripts in \code{islandid} islands.
      }
      
      \item{\code{signature(islandid = "missing", txid = "character", genomeDB = "annotatedGenome")}}{
	%%  ~~describe this method here~~
	Retrieve lengths from \code{genomeDB} for \code{txid} transcripts.
      }
      
      \item{\code{signature(islandid = "missing", txid = "missing", genomeDB = "annotatedGenome")}}{
	%%  ~~describe this method here~~
	Retrieve or calculate lengths for all transcripts in the annotated
	genome \code{genomeDB}.
      }
    }
  }
  \details{
    When called for the first time lengths are calculated and stored in
    the object \code{genomeDB}. Subsequent calls refer to these computed
    values.
  }
  \value{
   Named numeric vector with transcript lengths. 
  }
  
\keyword{ ~~ other possible keyword(s) ~~ }
