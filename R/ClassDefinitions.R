require(methods)


setClass("simulatedSamples", representation("list"))

valid_simulatedSamples <- function(object) {
  msg <- NULL
  if (!all(sapply(object, function(z) 'simTruth' %in% names(z)))) msg <- "All list elements must have an element 'simTruth'"
  if (!all(sapply(object, function(z) 'simExpr' %in% names(z)))) msg <- "All list elements must have an element 'simExpr'"
  if (!all(sapply(object, function(z) class(z$simTruth)=='data.frame'))) msg <- "Element simExpr must be of class data.frame"
  if (!all(sapply(object, function(z) class(z$simExpr)=='ExpressionSet'))) msg <- "Element simExpr must be of class ExpressionSet"
  if(!(is.null(msg))) { TRUE } else { msg }
}

setValidity("simulatedSamples", valid_simulatedSamples)

setMethod("show", signature(object="simulatedSamples"), function(object) {
  cat("simulatedSamples object with ",length(object)," simulated datasets (",ncol(object[[1]]$simExpr)," samples each)\n",sep="")
  cat("- 'coef' gets true differences between group means (returns matrix)\n")
  cat("- 'exprs' gets estimated expressions (returns list of ExpressionSets)\n")
  cat("- 'mergeBatches' combines exprs with a given ExpressionSet (returns list of ExpressionSets)\n")
}
)

setMethod("coef", signature(object="simulatedSamples"), function(object) {
  ans <- do.call(cbind,lapply(object, function(z) z$simTruth[,'mean.1'] - z$simTruth[,'mean.2']))
  rownames(ans) <- rownames(object[[1]]$simTruth)
  ans
}
)


setMethod("exprs", signature(object="simulatedSamples"), function(object) {
  lapply(object, '[[', "simExpr")
}
)

setMethod("[", signature(x="simulatedSamples"), function(x, i, ...) {
  new("simulatedSamples", lapply(x,'[')[i])
}
)

setMethod("[", signature(x="simulatedSamples",i="missing",j="index",drop="missing"), function(x, i, j, ..., drop) {
  new("simulatedSamples", lapply(x, function(z) { z$simExpr <- z$simExpr[,j]; return(z) }))
}
)


###############################################################

setClass("procBam", representation(pbam = "GRanges", junx="GRanges", stranded = "logical", plus="GRanges", minus="GRanges", pjunx="GRanges", mjunx="GRanges"))

setMethod("show", signature(object="procBam"), function(object) {
  cat("procBam object created from",ifelse(object@stranded,"stranded","non-stranded"),"reads\n")
  if(object@stranded) {
    cat("Contains",length(object@plus),"ranges in the positive strand corresponding to",length(unique(values(object@plus)$names)),"unique read pairs\n")
    cat("and",length(object@minus),"ranges in the negative strand corresponding to",length(unique(values(object@minus)$names)),"unique read pairs\n")
  } else cat("Contains",length(object@pbam),"ranges corresponding to",length(unique(values(object@pbam)$names)),"unique read pairs\n")
}
          )


setMethod("getReads", signature(x="procBam"), function(x) x@pbam)


###############################################################

setClass("readDistrs", representation(lenDis = "array", stDis = "function"))

setMethod("show", signature(object="readDistrs"), function(object) {
  cat("readDistrs object\n\n")
  cat("Insert size distribution (only first few shown)\n")
  show(head(object@lenDis,n=6))
  cat("...\n")
  cat("Read start cumulative distribution function in slot stDis\n")
}
)

setMethod("plot", signature(x="readDistrs"), function(x, y, ...) {
  if (missing(y)) stop("Specify type of plot: 'fragLength' or 'readSt'")
  args <- list(...)
  if ('col' %in% names(args)) col <- args$col else col <- 1
  if ('lty' %in% names(args)) lty <- args$lty else lty <- 1
  if ('lwd' %in% names(args)) lwd <- args$lwd else lwd <- 1
  if (y=='fragLength') {
    n <- as.numeric(names(x@lenDis))
    x2plot <- double(max(n)-min(n)+1); names(x2plot) <- min(n):max(n)
    x2plot[names(x@lenDis)] <- x@lenDis
    x <- as.numeric(names(x2plot))
    if ('xlim' %in% names(args)) xlim <- args$xlim else xlim <- range(x)
    y2plot <- x2plot/sum(x2plot)
    if ('ylim' %in% names(args)) ylim <- args$ylim else ylim <- c(0,max(y2plot))
    plot(x,y2plot,type='l',xlab='Fragment length',ylab='Proportion of reads',xlim=xlim,ylim=ylim,col=col,lty=lty,lwd=lwd)
  } else if (y=='readSt') {
    s <- seq(0,1,by=0.02)
    probs <- diff(x@stDis(s))
    if ('ylim' %in% names(args)) ylim <- args$ylim else ylim <- c(0,max(probs[!is.na(probs)]))
    plot(NA,NA,xlim=c(0,1),ylim=ylim,xlab='Read start (relative to transcript length)',ylab='Density')
    segments(s[-length(s)],probs,s[-1],col=col,lty=lty,lwd=lwd)
    segments(s,c(0,probs),s,c(probs,0),col=col,lty=lty,lwd=lwd)
  } else {
    stop("Second argument must be either 'fragLength' or 'readSt'")
  }
}
)


setMethod("lines", signature(x="readDistrs"), function(x, ...) {
  args <- list(...)
  y <- args[[1]]
  if ('col' %in% names(args)) col <- args$col else col <- 1
  if ('lty' %in% names(args)) lty <- args$lty else lty <- 1
  if ('lwd' %in% names(args)) lwd <- args$lwd else lwd <- 1
  if (y=='fragLength') {
    n <- as.numeric(names(x@lenDis))
    x2plot <- double(max(n)-min(n)+1); names(x2plot) <- min(n):max(n)
    x2plot[names(x@lenDis)] <- x@lenDis
    lines(as.numeric(names(x2plot)),x2plot/sum(x2plot),col=col,lty=lty,lwd=lwd)
  } else if (y=='readSt') {
    s <- seq(0,1,by=0.02)
    probs <- diff(x@stDis(s))
    segments(s[-length(s)],probs,s[-1],col=col,lty=lty,lwd=lwd)
    segments(s,c(0,probs),s,c(probs,0),col=col,lty=lty,lwd=lwd)
  } else {
    stop("Second argument must be either 'fragLength' or 'readSt'")
  }
}
)



###############################################################


setClass("pathCounts", representation(counts = "list", denovo = "logical", stranded="logical"))

valid_pathCounts <- function(object) {
  msg <- NULL
#validity checks go here
  if(object@stranded & length(object@counts)<2) msg <- "Stranded pathCounts must contain plus and minus counts"
  if(!(is.null(msg))) { TRUE } else { msg }
}

setValidity("pathCounts", valid_pathCounts)
setMethod("show", signature(object="pathCounts"), function(object) {
  if(object@denovo) {
    if(object@stranded) {
      cat("Stranded denovo pathCounts object with",length(object@counts[[1]])+length(object@counts[[2]]),"islands and", sum(!unlist(lapply(object@counts[['plus']], is.null))),"non zero positive islands and ", sum(!unlist(lapply(object@counts[['minus']], is.null))),"non zero negative islands\n")
          } else cat("Non-stranded denovo pathCounts object with",length(object@counts[[1]]),"islands and",sum(!unlist(lapply(object@counts[[1]], is.null))),"non zero islands.\n") 
  } else {
    if(object@stranded) {
      cat("Stranded known pathCounts object with",length(object@counts[[1]])+length(object@counts[[2]]),"islands and", sum(!unlist(lapply(object@counts[['plus']], is.null))),"non zero positive islands and ", sum(!unlist(lapply(object@counts[['minus']], is.null))),"non zero negative islands\n")
    } else cat("Non-stranded known pathCounts object with",length(object@counts[[1]]),"islands and", sum(!unlist(lapply(object@counts[[1]], is.null))),"non zero islands.\n")
  }
})



###############################################################


setClass("annotatedGenome", representation(islands = "GRangesList", transcripts = "list", knownVars = "list", exon2island = "data.frame", exonsNI="GRanges",  aliases="data.frame", genomeVersion="character", dateCreated="Date", txLength="numeric", denovo="logical"))
valid_annotatedGenome <- function(object) {
  msg <- NULL
  if (!(all(c('seqnames','start','end','width','island') %in% names(object@exon2island)))) msg <- "Incorrect column names in 'exon2island'"
  if(!(all(c('tx_id','tx_name','gene_id','tx','island_id') %in% names(object@aliases)))) msg <- "Incorrect column names in 'aliases'"
  if (is.null(msg)) { TRUE } else { msg }
}


setValidity("annotatedGenome", valid_annotatedGenome)
setMethod("show", signature(object="annotatedGenome"), function(object) {
  if(object@denovo) {
    cat("Denovo annotatedGenome object with",length(object@islands),"gene islands,", length(unique(unlist(lapply(object@transcripts, names)))),"transcripts and",nrow(object@exon2island),"exons.\n")
  } else cat("Known annotatedGenome object with",length(object@islands),"gene islands,", length(unique(unlist(lapply(object@transcripts, names)))),"transcripts and",nrow(object@exon2island),"exons.\n")
  cat("Genome version:",object@genomeVersion,"\n")
  cat("Date created:", as.character(object@dateCreated),"\n")
}
          )
