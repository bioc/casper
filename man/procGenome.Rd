\name{procGenome}
\alias{procGenome}
\alias{procGenome,TranscriptDb-method}
\alias{procGenome,GRanges-method}
\alias{procGenome,TranscriptDb,ANY-method}
\alias{procGenome,GRanges,ANY-method}
\alias{createDenovoGenome}
\title{
  Create an annotatedGenome object that stores information about genes and transcripts
}
\description{
  \code{procGenome} processes annotations for a given transcriptome,
  either from a \code{TranscriptDb} object created by \code{GenomicFeatures} package
  (e.g. from UCSC) or from a user-provided \code{GRanges} object (e.g. by
  importing a gtf file output by Cufflinks RABT module with
  \code{import} from \code{rtracklayer} package).
  
  \code{createDenovoGenome} creates a de novo annotated genome by
  combining UCSC annotations and observed RNA-seq data.
}
\usage{
procGenome(genDB, genome, mc.cores = 1)

createDenovoGenome(reads, DB, minLinks=2,
maxLinkDist=1e+05, maxDist=1000, minConn=2, minJunx=3, minLen=12, mc.cores=1)
}
\arguments{
  \item{genDB}{Either annotations downloaded from UCSC
    (e.g. genDB<-makeTranscriptDbFromUCSC(genome="hg19",
    tablename="refGene")) or \code{GRanges} object (e.g. genDB <-
    import('transcripts.gtf')). See details.}
  \item{genome}{Character indicating genome version (e.g. "hg19",
    "dm3")}
  \item{mc.cores}{Number of cores to use in parallel processing
    (multicore package required)}
  \item{DB}{\code{annotatedGenome} object, as returned by
    \code{procGenome}}
  \item{minLinks}{Minimum number of reads joining two
    exons to merge their corresponding genes}
  \item{maxLinkDist}{Maximum distance between two exons to merge
    their correspondin genes. A value of \code{0} disables this option.}
  \item{maxDist}{Maximum distance between two exons with reads joining
    them to merge their corresponding genes.}
  \item{minConn}{Minimum number of fragments connecting a new exon to an
    annotated one to add to denovo genome.}
  \item{minJunx}{Minimum number of junctions needed to redefine an annotated
    exon's end or start.}
  \item{minLen}{Minimum length of a junction to consider as a putative intron.}
  \item{reads}{Processed reads stored in a \code{RangedData}, as
    returned by \code{procBam}}
}
\section{Methods}{
\describe{
  \item{\code{signature(genDB = "transcriptDb")}}{
    \code{genDB} is usually obtained with a call to
    \code{makeTranscriptDbFromUCSC} (package \code{GenomicFeatures}),
    e.g. genDB<-makeTranscriptDbFromUCSC(genome="hg19", tablename="refGene")
  }
  \item{\code{signature(genDB = "GRanges")}}{
    \code{genDB} stores information about all transcripts and their
    respective exons. Chromosome, start, end and strand are stored as
    usual in GRanges objects. \code{genDB} must have a column named
    \code{"type"} taking the value \code{"transcript"} for rows
    corresponding to transcript and \code{"exon"} for rows corresponding
    to exons. It must also store transcript and gene ids. For instance, Cufflinks RABT
    module creates a gtf file with information formatted in this manner
    for known and de novo predicted isoforms.
  }
}}

\details{
  These functions create the annotation objects that are needed for subsequent
  functions. Typically these objects are created only once for a set of
  samples.

  If interested in quantifying expression for known transcripts
  only, one would typically use \code{procGenome} with a
  \code{TranscriptDb} from the usual Bioconductor
  annotations.
  It is also possible to provide \code{procGenome} with any desired set
  of transcripts, so that one can use other genome annotations or even
  de novo transcripts found by some third-party software. These
  annotations should be provided as a \code{GRanges} object. For
  instance, one could use the transcripts found by the Cufflinks RABT
  module. These transcripts are saved in gtf files, which can be
  imported into Bioconductor as \code{GRanges} with function \code{import} (package \code{rtracklayer}).
  \code{procGenome} is intended to be used in combination with
  \code{calcExp}, which quantifies expression for a set of given transcripts.
  
  \code{createDenovoGenome} performs a de novo transcript
  discovery following our own strategy, and is intended to be used in
  combination with \code{calcDenovo}, which performs a model search
  additional to quantifying expression. It is convenient to create a
  single genome by combining information from all samples, as this
  increases the power to detect new exons and
  guarantees that the same exons and gene islands are used across samples.
}
\value{
  Object of class \code{annotatedGenome}.
}
\author{
  Camille Stephan-Otto Attolini
}
\seealso{
  See \code{annotatedGenome-class} for a description of the class.
  To obtain the island id for an entrez gene id or a transcript id see help(getIsland).
  To extract the transcripts for a entrez or island id see
  help(transcripts).
  See \code{splitGenomeByLength} for splitting an \code{annotatedGenome}
  according to gene length.
}
\examples{
##Known transcripts from UCSC
## genDB<-makeTranscriptDbFromUCSC(genome="hg19", tablename="refGene")
## hg19DB <- procGenome(genDB, "hg19")

##Transcripts inferred by Cufflinks-RABT module
##genDB.Cuff <- import('transcripts.gtf')
##hg19DB.Cuff <- procGenome(genDB.Cuff, genome='hg19')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~annotation }
