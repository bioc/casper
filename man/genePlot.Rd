\name{genePlot}
\alias{genePlot}
\alias{genePlot-methods}
\alias{genePlot,CompressedIRangesList,ANY,ANY,ANY,ANY-method}
\alias{genePlot,IRanges,ANY,ANY,ANY,ANY-method}
\alias{genePlot,GRangesList,ANY,ANY,ANY,ANY-method}
\alias{genePlot,IRangesList,ANY,ANY,ANY,ANY-method}
\alias{genePlot,GRanges,ANY,ANY,ANY,ANY-method}
\alias{genePlot,missing,character,annotatedGenome,GRanges,ExpressionSet-method}
\alias{genePlot,missing,character,annotatedGenome,missing,missing-method}
\alias{genePlot,missing,character,annotatedGenome,procBam,ExpressionSet-method}
\alias{genePlot,missing,character,annotatedGenome,procBam,missing-method}
\title{
  Plot exon structure for each transcript of a given gene.
}
\description{
  Plot exon structure for each transcript of a given gene. Optionally, aligned reads can be
  added to the plot.
}
\usage{
genePlot(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='',
ylab='', xlim, cex=1, yaxt='n', col, ...)
}
\arguments{
  \item{generanges}{Object containing the ranges with start/end of each exon.}
  \item{islandid}{If \code{generanges} is not specified, transcripts are
  obtained from island \code{islandid} from the annotated genome \code{genomeDB}.}
  \item{genomeDB}{Annotated genome produced with the "procGenome" function}
  \item{reads}{\code{pbam} object with aligned reads. This is an
    optional argument.}
  \item{exp}{\code{ExpressionSet} object with expression values, as
    returned by \code{calcExp}. This is an optional argument.}
  \item{names.arg}{Optionally, indicate the names of each transcript.}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{xlim}{x-axis limits, defaults to start of 1st exon and end of
    last exon}
  \item{cex}{Character expansion}
  \item{yaxt}{The y-axis in the plot has no interpretation, hence by
    default it is not displayed.}
  \item{col}{Either single color or vector of colors to be used to draw
    each transcript. Defaults to rainbow colors.}
  \item{...}{Other arguments to be passed on to \code{plot}.}
}
\section{Methods}{
\describe{
  \item{\code{signature(generanges="CompressedIRangesList", islandid="ANY", genomeDB="ANY", reads="ANY", exp="ANY")}}{
    Plots a set of transcripts. Each element in the \code{generanges} corresponds to a
    transcript. Each transcript should contain exon start/end positions.
  }
  \item{\code{signature(generanges="IRanges", islandid="ANY", genomeDB="ANY", reads="ANY", exp="ANY")}}{
    Plots a single transcript. Each range indicates the start/end of a single exon.
  }
  \item{\code{signature(generanges="IRangesList", islandid="ANY", genomeDB="ANY", reads="ANY", exp="ANY")}}{
    Plots a set of transcripts. Each element in the \code{generanges} corresponds to a
    transcript. Each transcript should contain exon start/end positions.
  }
   \item{\code{signature(generanges="GRangesList", islandid="ANY",
       genomeDB="ANY", reads="ANY", exp="ANY")}}{
         Plots a set of transcripts. Each element in the
	 \code{generanges} corresponds to a
	 transcript. Each transcript should contain exon start/end
	 positions.
	   }
  \item{\code{signature(generanges="GRanges", islandid="ANY", genomeDB="ANY", reads="ANY", exp="ANY")}}{
    Plots a set of transcripts. Each space in \code{generanges} corresponds to a
    transcript. Each transcript should contain exon start/end positions.
  }
  \item{\code{signature(generanges="missing", islandid="character", genomeDB="annotatedGenome", reads="GRanges", exp="ExpressionSet")}}{
    Plots all transcripts stored in \code{genomeDB} for island with
    identifier \code{islandid}. Individual reads are added to the plot
    (\code{reads} contains start/end of individual read fragments).
  }
  \item{\code{signature(generanges="missing", islandid="character", genomeDB="annotatedGenome", reads="missing", exp="missing")}}{
    Plots all transcripts stored in \code{genomeDB} for island with identifier \code{islandid}.
  }
  \item{\code{signature(generanges="missing", islandid="character", genomeDB="annotatedGenome", reads="procBam", exp="missing")}}{
    Plots all transcripts stored in \code{genomeDB} for island with
    identifier \code{islandid}. Individual reads are added to the plot
    (\code{reads} contains start/end of individual read fragments).
  }
  \item{\code{signature(generanges="missing", islandid="character", genomeDB="annotatedGenome", reads="procBam", exp="ExpressionSet")}}{
    Plots all transcripts stored in \code{genomeDB} for island with
    identifier \code{islandid}. Individual reads and estimated
    expression are added to the plot
    (\code{reads} contains start/end of individual read fragments).
  }
}}
\value{
  A plot is produced.
}
\author{
  Camille Stephan-Otto Attolini, David Rossell
}
\examples{
data(hg19DB)

#Plot an IRangesList
txs <- transcripts(entrezid="27",genomeDB=hg19DB)
genePlot(txs)

#Equivalently, indicate islandid
islandid <- getIsland(entrezid="27",genomeDB=hg19DB)
genePlot(islandid=islandid, genomeDB=hg19DB)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
