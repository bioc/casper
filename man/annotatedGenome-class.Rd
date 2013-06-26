\name{annotatedGenome-class}
\Rdversion{1.1}
\docType{class}
\alias{annotatedGenome}
\alias{annotatedGenome-class}
\alias{show,annotatedGenome-method}

\title{Class "annotatedGenome"}
\description{
  Stores an annotated genome, either known or build de novo by combining
  a known genome with observed RNA-seq data.
}
\section{Objects from the Class}{
  Objects are typically created with a call to
  \code{procGenome} (for known genomes) or to
  \code{createDenovoGenome} (for de novo genomes).
}
\section{Slots}{
  \describe{
    \item{islands}{\code{GRangesList} object with elements corresponding
      to gene islands. It indicates the
      start/end/name of each exon contained in the island}
    \item{transcripts}{Each element in the list corresponds to a gene island. It indicates the exons contained in each known variant.}
    \item{exon2island}{\code{data.frame} indicating the chromosome,
      start and end of each exon, and its corresponding gene island.}
    \item{exonsNI}{\code{GRanges} indicating the chromosome,
      start/end and id of each exon}
    \item{aliases}{\code{data.frame} indicating the aliases for each
      known transcript, i.e. transcripts having the exact same sequence
      of exons (as indicated in the UCSC database).}
    \item{genomeVersion}{Character indicating the genome version from
      which the object was build, e.g. "hg19"}
    \item{dateCreated}{Character indicating the date when the object was
      created. UCSC genomes chance from time to time, so that an "hg19"
      genome from Jan 2012 may not be exactly the same as in Dec 2012. }
    \item{denovo}{Logical variable. \code{FALSE} indicates that the
      object was created using available annotation only. \code{TRUE}
      indicates that new exons/islands were added based on the data
      observed in a particular RNA-seq experiment.}
    \item{txLength}{Numeric vector storing transcript lengths.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "annotatedGenome")}: Displays general information about the object. }
	 }
}
\author{
Camille Stephan-Otto Attolini
}
\seealso{
  \code{\link{procGenome}} and \code{\link{createDenovoGenome}} to
  create \code{annotatedGenome} objects.
}
\examples{
showClass("annotatedGenome")
}
\keyword{classes}
