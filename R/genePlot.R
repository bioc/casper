buildGene<-function(genomeDB, islandid, txs){
    rd<-lapply(txs, function(x) genomeDB@islands[[islandid]][names(genomeDB@islands[[islandid]]) %in% names(x),])
    rdl<-GRangesList(rd)
    rdl
}

setMethod("genePlot", signature(generanges='missing',islandid='character',genomeDB='annotatedGenome',reads='missing',exp='missing') , function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
  txs <- transcripts(islandid=islandid,genomeDB=genomeDB)
  genePlot(txs, names.arg=names.arg, xlab=xlab, ylab=ylab, xlim=xlim, cex=cex, yaxt=yaxt, col=col, ...)
}
)

setMethod("genePlot", signature(generanges='missing',islandid='character',genomeDB='annotatedGenome',reads='procBam',exp='ExpressionSet') , function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
  txs <- transcripts(islandid=islandid, genomeDB=genomeDB)
  txexp<-exprs(exp)[names(txs),]
  #if(sum(txexp)>0) txexp<-(txexp/sum(txexp))*100
  gene<-buildGene(txs=txs, genomeDB=genomeDB, islandid=islandid)
  names(gene) <- paste(names(gene),' (expr=',round(txexp[names(gene)],3),')',sep='')
  chr <- getChr(islandid=islandid, genomeDB=genomeDB)
  reads@pbam <- reads@pbam[seqnames(reads@pbam) %in% chr]
  rangesPlot(x=reads, gene, exonProfile=FALSE)
}
)

setMethod("genePlot", signature(generanges='missing',islandid='character',genomeDB='annotatedGenome',reads='GRanges',exp='ExpressionSet') , function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
  txs <- transcripts(islandid=islandid, genomeDB=genomeDB)
  txexp<-exprs(exp)[names(txs),]
  #if(sum(txexp)>0) txexp<-(txexp/sum(txexp))*100
  gene<-buildGene(txs=txs, genomeDB=genomeDB, islandid=islandid)
  rangesPlot(x=reads[as.character(genomeDB@exon2island[rownames(genomeDB@exon2island)==names(gene[[1]])[1],]$space)], gene, exonProfile=FALSE)
  #txexp
}
)

#setMethod("plotExpr",signature(gene="denovoGeneExpr"),
#  function(gene, minProbExpr, minExpr, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
#    n <- names(variants(gene))
#    n[nchar(n)>30] <- paste("Variant ",1:sum(nchar(n)>30),sep="")
#    names(n) <- names(variants(gene))
#    v <- variantMargExpr(gene,minProbExpr=minProbExpr,minExpr=minExpr)
#    names.arg <- n[rownames(v)]
#    names.arg <- paste(names.arg,' (Rel. expr.=',round(v[,'expr'],3),'; P(expressed)=',round(v[,'probExpressed'],3),')',sep='')
#    names(names.arg) <- n[rownames(v)]
#    gene <- variants(gene)[rownames(v)]
#    o <- names(n)[names(n) %in% names(gene)]
#    gene <- gene[o]
#    names.arg <- names.arg[n[o]]
#    start(gene) <- round(start(gene)/1000)
#    end(gene) <- round(end(gene)/1000)
#    if (missing(xlim)) xlim <- c(min(start(gene)),max(end(gene)))
#    genePlot(gene, names.arg=names.arg, xlab=xlab, ylab=ylab, xlim=xlim, cex=1, yaxt='n', col, ...)
#  }
#)


setMethod("genePlot",signature(generanges='IRanges'),
  function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
    if (missing(xlim)) xlim <- range(c(start(generanges),end(generanges)))
    if (missing(col)) col <- 1
    plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab=xlab,ylab=ylab,yaxt=yaxt,...)
    nsplice <- 1
    y0 <- seq(1/(1+nsplice),1-1/(1+nsplice),by=1/(1+nsplice))
    x0 <- start(generanges)
    x1 <- end(generanges)
    if (length(x0)>1) segments(x0=x1[-length(x1)],x1=x0[-1],y0=y0,y1=y0)
    segments(x0=x0,x1=x1,y0=y0+.33/(1+nsplice),y1=y0+.33/(1+nsplice),col=col)
    segments(x0=x0,x1=x1,y0=y0-.33/(1+nsplice),y1=y0-.33/(1+nsplice),col=col)
    segments(x0=x0,x1=x0,y0=y0-.33/(1+nsplice),y1=y0+.33/(1+nsplice),col=col)
    segments(x0=x1,x1=x1,y0=y0-.33/(1+nsplice),y1=y0+.33/(1+nsplice),col=col)
  }
)


setMethod("genePlot",signature(generanges='GRanges'),
  function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
    genePlot(ranges(generanges), names.arg=names.arg, xlab=xlab, ylab=ylab, xlim=xlim, cex=cex, yaxt=yaxt, col=col, ...)
  }
)

setMethod("genePlot",signature(generanges='IRangesList'),
  function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
    if (missing(xlim)) xlim <- range(unlist(c(start(generanges),end(generanges))))
    if (missing(names.arg)) names.arg <- names(generanges)
    plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab=xlab,ylab=ylab,yaxt=yaxt,...)
    nsplice <- length(generanges)
    if (missing(col)) col<-rep(rainbow(min(nsplice, 10)), ceiling(nsplice/10)) else if (length(col)!=nsplice) col <- rep(col[1],nsplice)
    y0 <- seq(1/(1+nsplice),1-1/(1+nsplice),by=1/(1+nsplice))
    for (i in nsplice:1) {
      x0 <- start(generanges[[i]])
      x1 <- end(generanges[[i]])
      if (length(x0)>1) segments(x0=x1[-length(x1)],x1=x0[-1],y0=y0[i],y1=y0[i], col=col[i])
      segments(x0=x0,x1=x1,y0=y0[i]+.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
      segments(x0=x0,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]-.33/(1+nsplice), col=col[i])
      segments(x0=x0,x1=x0,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
      segments(x0=x1,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
    }
    if (!is.null(names.arg)) text(rep(xlim[1],length(y0)),y0+.5*(y0[2]-y0[1]),names.arg,cex=cex,adj=0)
  }
)

setMethod("genePlot",signature(generanges='GRangesList'),
  function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
    if (missing(xlim)) xlim <- range(c(start(unlist(generanges)),end(unlist(generanges))))
    if (missing(names.arg)) names.arg <- names(generanges)
    plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab=xlab,ylab=ylab,yaxt=yaxt,...)
    nsplice <- length(generanges)
    if (missing(col)) col<-rep(rainbow(min(nsplice, 10)), ceiling(nsplice/10)) else if (length(col)!=nsplice) col <- rep(col[1],nsplice)
    y0 <- rev(seq(1/(1+nsplice),1-1/(1+nsplice),by=1/(1+nsplice)))
    for (i in 1:nsplice) {
      x0 <- start(generanges[[i]])
      x1 <- end(generanges[[i]])
      if (length(x0)>1) segments(x0=x1[-length(x1)],x1=x0[-1],y0=y0[i],y1=y0[i], col=col[i])
      segments(x0=x0,x1=x1,y0=y0[i]+.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
      segments(x0=x0,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]-.33/(1+nsplice), col=col[i])
      segments(x0=x0,x1=x0,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
      segments(x0=x1,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
    }
    if (!is.null(names.arg)) text(rep(xlim[1],length(y0)),y0+.5*(y0[1]-y0[2]),names.arg,cex=cex,adj=0)
  }
)




#setMethod("rangesPlot",signature(x='IRanges',gene='IRanges'),
#  function(x, gene, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
#    txStart <- min(start(gene))
#    txEnd <- max(end(gene))
#    if (missing(xlim)) xlim <- c(txStart,txEnd)
#    x <- x[start(x)>=xlim[1] & end(x)<=xlim[2]]
#    if (length(x)==0) {
#      warning('There are no reads in the specified region')
#    } else {
#      x <- x[order(start(x))]
#      cover <- coverage(x)
#      cover <- seqselect(cover,(cover@lengths[1]+1):length(cover))
#      layout(c(2,1),heights=heights)
#      par(oma=c(1,0,0,0),mar=c(2.1,.1,.1,.1))
#      genePlot(gene=gene, xlab='', ylab='', xlim=xlim, ...)
#      sel <- width(x)<=maxFragLength
#      if (exonProfile) {
#        st <- c(start(x[sel]),start(x[!sel]),end(x[!sel])); en <- c(end(x[sel]),start(x[!sel]),end(x[!sel]))
#        pr <- coverage(IRanges(start=st,end=en))
#        pr <- seqselect(pr, (pr@lengths[1]+1):length(pr))
#        minst <- min(st)
#        lines(minst:(minst+length(pr)-1),pr/max(pr),col='gray')
#      }
#      #Height for large reads (i.e. spanning over several exons)
#      plot(NA,NA,xlim=xlim,ylim=0:1,xlab='',ylab='',xaxt='n',yaxt='n',mgp=c(0,0,0))
#      x0 <- start(x); x1 <- end(x)
#      if (any(!sel)) {
#        y0 <- runif(sum(!sel))
#        segments(x0=x0[!sel],x1=x1[!sel],y0=y0,y1=y0)
#      }
#      #Height for short reads
#      if (any(sel)) {
#        y0 <- as.integer(cover)[start(x)[sel]-start(x)[1]+1] + runif(sum(sel))
#        y0 <- y0/max(y0)
#        segments(x0=x0[sel],x1=x1[sel],y0=y0,y1=y0,col='gray')
#      }
#    }
#  }
#)


setMethod("rangesPlot",signature(x='IRanges'),
  function(x, gene, exonProfile=TRUE, maxFragLength=500, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    x <- x[start(x)>=xlim[1] & end(x)<=xlim[2]]
    if (length(x)==0) {
      warning('There are no reads in the specified region')
    } else {
      x <- x[order(start(x))]
      cover <- coverage(x)
      cover <- seqselect(cover,(cover@lengths[1]+1):length(cover))
      layout(c(2,1),heights=heights)
      par(oma=c(1,0,0,0),mar=c(2.1,.1,.1,.1))
      genePlot(gene=gene, xlab='', ylab='', xlim=xlim, ...)
      sel <- width(x)<=maxFragLength
      if (exonProfile) {
        st <- c(start(x[sel]),start(x[!sel]),end(x[!sel])); en <- c(end(x[sel]),start(x[!sel]),end(x[!sel]))
        pr <- coverage(IRanges(start=st,end=en))
        pr <- seqselect(pr, (pr@lengths[1]+1):length(pr))
        minst <- min(st)
        lines(minst:(minst+length(pr)-1),pr/max(pr),col='gray')
      }
      #Height for large reads (i.e. spanning over several exons)
      plot(NA,NA,xlim=xlim,ylim=0:1,xlab='',ylab='',xaxt='n',yaxt='n',mgp=c(0,0,0))
      x0 <- start(x); x1 <- end(x)
      if (any(!sel)) {
        y0 <- as.integer(cover)[start(x)[!sel]-min(start(x))+1] + runif(sum(!sel))
        y0 <- y0/max(y0)
       # y0 <- runif(sum(!sel))
        segments(x0=x0[!sel],x1=x1[!sel],y0=y0,y1=y0)
      }
      #Height for short reads
      if (any(sel)) {
        y0 <- as.integer(cover)[start(x)[sel]-min(start(x))+1] + runif(sum(sel))
        y0 <- y0/max(y0)
        segments(x0=x0[sel],x1=x1[sel],y0=y0,y1=y0,col='gray')
      }
    }
  }
)


setMethod("rangesPlot",signature(x='GRanges',gene='IRanges'),
  function(x, gene, exonProfile=TRUE, maxFragLength=500, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    if (missing(xlim)) xlim <- c(min(start(gene)),max(end(gene)))
    rangesPlot(unlist(ranges(x)),gene=gene,exonProfile=exonProfile,maxFragLength=maxFragLength,xlab=xlab,ylab=ylab,xlim=xlim,heights=heights,...) 
  }
)

setMethod("rangesPlot",signature(x='GRanges',gene='GRangesList'),
function(x, gene, exonProfile=TRUE, maxFragLength=500, xlab='', ylab='', xlim, heights=c(2,1), ...) {
  if (missing(xlim)) xlim <- c(min(unlist(start(gene))),max(unlist(end(gene))))
  rangesPlot(unlist(ranges(x)),gene=gene,exonProfile=exonProfile,maxFragLength=maxFragLength,xlab=xlab,ylab=ylab,xlim=xlim,heights=heights,...)  
}
)



setMethod("rangesPlot",signature(x='procBam'),
  function(x, gene, exonProfile=TRUE, maxFragLength=500, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    if (missing(xlim)) xlim <- range(unlist(start(gene)))
    x <- getReads(x)
    x <- x[start(x)>=xlim[1] & end(x)<=xlim[2],]
    if(sum(duplicated(names(x)))==0) names(x) <- sub("\\..*", "", names(x))
    if (length(x)==0) {
      warning('There are no reads in the specified region')
    } else {
      layout(c(2,1),heights=heights)
      par(oma=c(1,0,0,0),mar=c(2.1,.1,.1,.1))
      genePlot(gene=gene, xlab='', ylab='', xlim=xlim, ...)
      tab <- table(names(x))
      sel <- names(x) %in% c(names(tab[tab>2]), findLongInserts(x, minsize=maxFragLength))
      notsel <- !sel
      if (exonProfile) {
        st <- c(start(x)[sel],start(x)[!sel],end(x)[!sel]); en <- c(end(x)[sel],start(x)[!sel],end(x)[!sel])
        pr <- coverage(IRanges(start=st,end=en))
        pr <- seqselect(pr, (pr@lengths[1]+1):length(pr))
        minst <- min(st)
        lines(minst:(minst+length(pr)-1),pr/max(pr),col='gray')
      }
      #Height for large reads (i.e. spanning over several exons)
      plot(NA,NA,xlim=xlim,ylim=0:1,xlab='',ylab='',xaxt='n',yaxt='n',mgp=c(0,0,0))
      x0 <- start(x); x1 <- end(x)
      y0 <- as.numeric(names(x))/max(as.numeric(names(x)))
      if (any(notsel)) {
        segments(x0=x0[notsel],x1=x1[notsel],y0=y0[notsel])
      }
      #Long reads
      if (any(sel)) {
        longst <- tapply(x0[sel],INDEX=names(x)[sel],FUN=min)
        longend <- tapply(x1[sel],INDEX=names(x)[sel],FUN=max)
        longy0 <- tapply(y0[sel],INDEX=names(x)[sel],FUN=function(z) z[1])
        segments(x0=longst,x1=longend,y0=longy0,col=4,lty=3)
        segments(x0=x0[sel],x1=x1[sel],y0=y0[sel],col='red',lwd=3)
      }
    }
  }
)



findLongInserts <- function(x, minsize) {
  #Select long paired end inserts in a GRanges object ordered by pair id (e.g. as returned by procBam in package casper)
  # - x: GRanges, column id must indicate read pair id. Assumed to be ordered according to read id & fragment start
  # - minsize: minimum length for an insert to be considered long. Defaults to 99% percentile of observed insert widths
  # Returns: ids of long inserts
  
  w <- end(x)[-1] - start(x)[-length(x)]
  sel <- names(x)[-length(x)] == names(x)[-1]
  w[!sel] <- 0
  if (missing(minsize)) minsize <- quantile(w[sel],probs=.99)
  #Select long inserts
  idsel <- names(x)[which(w>minsize)]
  sel <- which(names(x) %in% idsel)
  chr <- seqnames(x)[sel]; st <- start(x)[sel]; en <- end(x)[sel]; id <- names(x)[sel]
  leftend <- which(c(TRUE,id[-1]!=id[-length(id)]))
  rightend <- c(leftend[-1]-1,length(id))
  st <- st[leftend]; en <- en[rightend]; chr <- chr[leftend]; id <- id[leftend]
  id[st<en]
}
