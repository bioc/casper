setGeneric("findNewExons", function(pbam, DB, cov, minConn, minJunx, minLen, mc.cores) standardGeneric("findNewExons"))
setMethod("findNewExons", signature(pbam='list'),
          function(pbam, DB, cov, minConn, minJunx, minLen, mc.cores) {
            if(mc.cores>1) require(parallel)
            if(mc.cores>1) {
              if ('parallel' %in% loadedNamespaces()) {
                nex <- parallel::mclapply(pbam, function(x) findNewExons(x, DB=DB, cov=cov, minConn=minConn, minJunx=minJunx, minLen=minLen), mc.cores=mc.cores)
              } else stop('parallel library has not been loaded!')
            } else {
              nex <- lapply(pbam, function(x) findNewExons(x, DB=DB, cov=cov, minConn=minConn, minJunx=minJunx, minLen=minLen))
            }
            nex <- casper:::mergeGRanges(nex)
            names(nex) <- 1:length(nex)
            nex
          }
          )
setMethod("findNewExons", signature(pbam='procBam') ,
          function(pbam, DB, cov, minConn, minJunx, minLen) {
            nex <- findNewExonsF(pbam=pbam, DB=DB, cov=cov, minConn=minConn, minJunx=minJunx, minLen=minLen)
            names(nex) <- 1:length(nex)
            nex
          }
          )
 
setGeneric("assignExons2Gene", function(reads, chrs, ...) standardGeneric("assignExons2Gene"))
setMethod("assignExons2Gene", signature(chrs='list'),
          function(reads, chrs, mc.cores, DB, ...){
            ans <- parallel::mclapply(chrs, function(x) assignExons2Gene(reads=reads, DB=DB, chrs=x, mc.cores=mc.cores, ...), mc.cores=mc.cores)
            misschr <- unlist(lapply(ans, function(x) as.character(unique(seqnames(x@exonsNI)))))
            allchr <- as.character(unique(seqnames(DB@exonsNI)))
            misschr <- allchr[!(allchr %in% misschr)]
            misschr <- lapply(misschr, function(x) formatChromo(subsetGenome(genomeDB=DB, chr=x)))
            ans <- mergeGenomeByChr(c(ans, misschr))
            ans
          })
setMethod("assignExons2Gene", signature(chrs='character', reads='procBam'),
          function(reads, chrs, DB, mc.cores, ...){
            assignExons2GeneF(chrs=chrs, reads=reads@pbam, DB=DB, ...)
          })
                   
formatChromo <- function(DB){
  chr <- as.character(unique(seqnames(DB@exonsNI)))
  isl <- DB@islands
  names(isl) <- paste(chr, names(isl), sep='.')
  txs <- DB@transcripts
  names(txs) <- names(isl)
  alia <- DB@aliases
  alia$island_id <- paste(chr, alia$island_id, sep=".")
  e2i <- DB@exon2island
  e2i$island <- paste(chr, e2i$island, sep='.')
  ans <- new("annotatedGenome", islands=isl, transcripts=txs, aliases=alia, exon2island=e2i, exonsNI=DB@exonsNI, denovo=T, dateCreated=DB@dateCreated, genomeVersion=DB@genomeVersion)
}

mergeGenomeByChr <- function(ans){
  isl <- lapply(ans, function(x) x@islands)
  isl <- do.call('c', isl)
  txs <- lapply(ans, function(x) x@transcripts)
  txs <- do.call('c', txs)
  alia <- lapply(ans, function(x) x@aliases[,c("tx_id","tx_name","gene_id","tx","island_id")])
  alia <- do.call('rbind', alia)
  e2i <- lapply(ans, function(x) x@exon2island)
  e2i <- do.call('rbind', e2i)
  exn <- lapply(ans, function(x) x@exonsNI)
  exn <- do.call('c', exn)
  new("annotatedGenome", islands=isl, transcripts=txs, aliases=alia, exon2island=e2i, exonsNI=exn, denovo=T, dateCreated=ans[[1]]@dateCreated, genomeVersion=ans[[1]]@genomeVersion)
}

findNewExonsF <- function(pbam, cov=NULL, DB, minConn=3, minJunx=3, minLen=12, mc.cores=1){
  junx <- pbam@junx
  pbam <- pbam@pbam
#  chrs <- as.character(unique(seqnames(pbam)@values))
#  DB <- subsetGenome(genomeDB=DB, chr=chrs)
  ## Find new putative exons by reads' islands
  if(missing(cov)) cov <- coverage(pbam)
  isl <- slice(cov, lower=1,  rangesOnly=T)
  islcov <- Views(cov, isl[names(cov)])

  cnt <- parallel::mclapply(names(isl), function(chr){
    viewApply(islcov[[chr]], sum)
  }, mc.cores=mc.cores)
  cnt <- unlist(cnt)
  len <- unlist(width(isl))
  cnt <- cnt/len
  isl <- GRanges(ranges=unlist(isl), seqnames=rep(names(isl), width(isl@partitioning)))
  isl <- isl[cnt>2]
  isls <- c(isl, DB@exonsNI)
  strand(isls) <- "*"
  isls <- reduce(isls)
  isls <- c(isls, DB@exonsNI)
  strand(isls) <- "*"
  isls <- disjoin(isls)
  
  ## Redefine known and new exons by junctions
  isls <- parallel::mclapply(unique(as.character(seqnames(DB@exonsNI)@values)), function(x){
    if(sum(seqnames(junx)==x)){
      y <- junx[seqnames(junx) == x]
      y <- y[width(y)>3]
      #z <- names(y) #paste(start(y), end(y), sep='.')
      zz <- values(y)$counts #table(z)
      names(zz) <- names(y)
      y <- y[names(zz)[zz>=minJunx]]
      values(y) <- NULL
      y <- disjoin(c(isls[seqnames(isls)==x], y))
      ans <- IRanges::unique(subsetByOverlaps(y, isls[seqnames(isls)==x]))
    } else ans <- disjoin(isls[seqnames(isls)==x])
    ans
  }, mc.cores=mc.cores)
  isls <- do.call('c', isls)
  names(isls) <- 1:length(isls)

  ## Find exons with no overlap with known exons
  kex <- subsetByOverlaps(isls, DB@exonsNI)
  nisl <- isls[!(isls %in% DB@exonsNI)]
  over <- findOverlaps(nisl, DB@exonsNI, maxgap=1)
  nex <- nisl[-queryHits(over)]

  ## Find exons to redefine (grow)
  lex <- nisl[!(names(nisl) %in% names(nex))]
  
  lex <- parallel::mclapply(unique(as.character(seqnames(DB@exonsNI)@values)), function(x){
    junx <- junx[seqnames(junx)==x]
    junx <- junx[width(junx)>5]
    scnt <- tapply(values(junx)$counts, start(junx), sum)
    ecnt <- tapply(values(junx)$counts, end(junx), sum)
    tmp <- lex[seqnames(lex)==x]
    cnts <- cbind(ecnt[as.character(start(tmp)-1)], scnt[as.character(end(tmp)+1)])
    cnts[is.na(cnts)] <- 0
    tmp <- tmp[cnts[,1]>=minJunx | cnts[,2]>=minJunx]
    tmp[width(tmp)>minLen]
  }, mc.cores=mc.cores)
  lex <- do.call('c', lex)
  ans <- c(kex, nex, lex)
  ans <- sort(ans)
  names(ans) <- 1:length(ans)
  return(ans)
}

assignExons2GeneF <- function(exons, DB, reads, chrs, maxDist=1000, minLinks=3, maxLinkDist=0, stranded=FALSE, mc.cores=1){
  #chrs <- as.character(unique(seqnames(reads)@values))
  reads <- subset(reads, as.character(seqnames(reads)) %in% chrs)
  DB <- subsetGenome(genomeDB=DB, chr=chrs)
  exons <- exons[seqnames(exons) %in% chrs]
  over <- findOverlaps(exons, DB@exonsNI)
  exid <- names(exons)[queryHits(over)]
  exkey <- split(exid, names(DB@exonsNI)[subjectHits(over)])
  exlen <- sapply(exkey, length)
  otxs <- unlist(DB@transcripts, recursive=F)
  names(otxs) <- sub("[0-9]+\\.", "", names(otxs))
  newTxs <- unlist(exkey[as.character(sprintf("%d",unlist(otxs)))])
  lenkey <- exlen[as.character(sprintf("%d",unlist(otxs)))]
  lenkey <- tapply(lenkey, rep(names(otxs), sapply(otxs, length)), sum)
  lenkey <- lenkey[names(otxs)]
  newTxs <- unname(as.integer(newTxs))
  newTxs <- split(newTxs, rep(names(lenkey), lenkey))
  newTxs <- newTxs[names(otxs)]
  newTxs <- split(newTxs, rep(names(DB@transcripts), sapply(DB@transcripts, length)))
  alljunx <- NULL
  exs <- exons[-queryHits(over)]
  if(length(exs)>0){
    #rea <- subsetByOverlaps(reads, exs)
  
#Select read ids in these exons
    #area <- reads[names(reads) %in% names(rea)]
    over <- findOverlaps(exons, reads)
    exs1 <- names(exons)[queryHits(over)]
    rea1 <- names(reads)[subjectHits(over)]
    rea1 <- sub("\\..*", "", rea1)
    len <- length(unique(rea1))
    if(length(unique(exs1))>1){
      ans <- .Call("joinExons", as.integer(exs1), as.integer(rea1), len)
      junx <- ans[[1]][ans[[2]]>=minLinks]
      junx <- strsplit(junx, split=".", fixed=T)
      names(junx) <- 1:length(junx)
      nalljunx <- rep(1:length(junx), unlist(lapply(junx, length)))
      alljunx <- unlist(junx)
      tmpex <- as.data.frame(exons)
      rownames(tmpex) <- names(exons)
      tmpex <- tmpex[alljunx,]
      alljunx <- as.integer(alljunx)
      names(alljunx) <- nalljunx
      pos <- (tmpex$end+tmpex$start)/2
      dis <- tapply(pos, names(alljunx), function(x) ifelse(length(x)>1, max(sort(x)[-length(x)]-sort(x)[-1]) < maxLinkDist, TRUE ))
      alljunx <- alljunx[names(alljunx) %in% as.character(names(dis)[dis])]
      if(sum(!(as.numeric(names(exs)) %in% unique(alljunx)))>0) {
        sing <- as.numeric(names(exs))[!(as.numeric(names(exs)) %in% unique(alljunx))]
        mname <- max(as.numeric(names(alljunx)))
        names(sing) <- (mname+1):(mname+length(sing)) 
        alljunx <- c(alljunx, sing)
      }
    } else {
      alljunx <- unique(exs1)
      names(alljunx) <- '1'
    }
  }

#Build gene islands
  oldexs <- unlist(newTxs, recursive=F)
  txids <- sub("[0-9]+\\.", "", names(oldexs))
  txids <- rep(txids, unlist(lapply(oldexs, length)))
  oldexs <- unlist(oldexs)
  names(oldexs) <- txids
  naj <- names(alljunx)
  alljunx <- as.integer(alljunx)
  names(alljunx) <- naj
  allexs <- c(oldexs, alljunx)
  islands <- casper:::makeIslands(allexs)

  nislands <- names(islands)
  islands <- paste(unique(as.character(seqnames(exons)@values)), islands, sep='.')
  names(islands) <- nislands
  values(exons)$island <- islands[names(exons)] 
  exon2gene <- as.data.frame(exons)
  rownames(exon2gene) <- names(exons)
  exon2gene$island <- islands[rownames(exon2gene)]
  islands <- GRanges(IRanges(exon2gene$start, exon2gene$end), seqnames=exon2gene$seqnames)
  names(islands) <- rownames(exon2gene)
  exon2island <- as.data.frame(exons)  
  exon2island$strand <- NULL

  oldstr <- as.character(strand(DB@islands@unlistData))
  oldstr[oldstr=='*'] <- "."
  sel <- match(islands, DB@islands@unlistData)
  isstr <- as.character(strand(islands))
  nstr <- as.character(oldstr)[sel]
  isstr[!is.na(sel)] <- nstr[!is.na(nstr)]
  isstr <- tapply(isstr, as.character(exon2gene$island), unique)
  mix <- sapply(isstr, function(x) sum(grepl(".", x, fixed=T))>0)
  mixx <- unlist(sapply(isstr, function(x) sum(grepl("-", x))>0 & sum(grepl("+", x, fixed=T))>0))
  mix <- mix | mixx
  plus <- sapply(isstr, function(x) sum(grepl("+", x, fixed=T))>0)
  minus <- sapply(isstr, function(x) sum(grepl("-", x, fixed=T))>0)
  fstr <- rep("*", length(plus))
  fstr[plus] <- "+"
  fstr[minus] <- "-"
  fstr[mix] <- "*"
  strand(islands) <- "+"
  islands <- split(islands, as.character(exon2gene$island[match(rownames(exon2gene), names(islands))]))
  
  strand(islands[fstr=='-']@unlistData) <- "-"
  strand(islands[fstr=='*']@unlistData) <- "*"
  islands[fstr=='-'] <- GRangesList(lapply(islands[fstr=='-'], rev))

  extxs <- unlist(newTxs, recursive=F)
  names(extxs) <- sub("[0-9]+\\.", "", names(extxs))
  sel <- unlist(lapply(extxs, "[", 1))
  sel <- match(as.character(sel), rownames(exon2island))
  tx2gene <- exon2gene$island[sel]
  transcripts <- vector(length=length(islands), mode="list")
  names(transcripts) <- names(islands)
  tmp <- split(extxs, tx2gene)
  transcripts[names(tmp)] <- tmp
  aliases <- DB@aliases
  aliases$island_id <- tx2gene[match(rownames(aliases),names(extxs))]
  aliases <- aliases[,!(colnames(aliases) %in% 'exid')]
  values(exons)$island <- NULL
  ans <- new("annotatedGenome", islands=islands, transcripts=transcripts, exon2island=exon2island, exonsNI=exons, aliases=aliases, genomeVersion=DB@genomeVersion, dateCreated=Sys.Date(), denovo=TRUE)
  ans <- subsetGenome(islands=names(ans@islands)[!(elementLengths(ans@islands)==1 & sapply(ans@transcripts, is.null))], genomeDB=ans)
  ans
}

mergeStrDenovo <- function(plus, minus){  
  nullplus <- sapply(plus@transcripts, is.null)
  if(sum(nullplus)>0) {
    newplus <- tapply(names(plus@islands[nullplus]@unlistData), rep(names(plus@islands)[nullplus], elementLengths(plus@islands[nullplus])), function(x) paste(x, collapse='.'))
  } else newplus <- NULL
  nullminus <- unlist(lapply(minus@transcripts, is.null))
  if(sum(nullminus)>0) {
    newminus <- tapply(names(minus@islands[nullminus]@unlistData), rep(names(minus@islands)[nullminus], elementLengths(minus@islands[nullminus])), function(x) paste(sort(x), collapse='.'))
  } else newminus <- NULL
  common <- NULL
  if(!is.null(newplus) & !is.null(newminus)) common <- names(newminus)[!(newminus %in% newplus)]
  allislands <- c(plus@islands, minus@islands[!(names(minus@islands) %in% common)])
  names(allislands) <- c(paste(names(plus@islands), "P", sep='.'), paste(names(minus@islands)[!(names(minus@islands) %in% common)], "M", sep="."))
  alltrans <- c(plus@transcripts, minus@transcripts[!(names(minus@islands) %in% common)])
  names(alltrans) <- c(paste(names(plus@islands), "P", sep='.'), paste(names(minus@islands)[!(names(minus@islands) %in% common)], "M", sep="."))
  exP <- plus@exonsNI
  values(exP)$island <- paste(values(exP)$island, ".P", sep="")
  exM <- minus@exonsNI
  values(exM)$island <- paste(values(exM)$island, ".M", sep="")
  allexonsNI <- unique(c(exP, exM))
  ex2is <- as.data.frame(allexonsNI, row.names=NULL)
  ex2is$id <- names(allexonsNI)
  alP <- plus@aliases
  alP$island_id <- paste(alP$island_id, ".P", sep="")
  alM <- minus@aliases
  alM$island_id <- paste(alM$island_id, ".M", sep="")
  alia <- rbind(alP, alM)
  ans <- new("annotatedGenome", aliases=alia, denovo=TRUE, exonsNI=allexonsNI, transcripts=alltrans, exon2island=ex2is, dateCreated=Sys.Date(), genomeVersion=plus@genomeVersion, islands=allislands)
  ans
}
createDenovoGenome <- function(reads, DB, minLinks=2,  maxLinkDist=1e+05, maxDist=1000, minConn=2, minJunx=3, minLen=12, mc.cores=1){}
#createDenovoGenome <- function(reads, DB, cov, minLinks=2, maxLinkDist=1e+05, maxDist=1000, minConn=2, minJunx=3, minLen=12, mc.cores=1){
#  cat("Finding new exons\n")
#  somex <- NULL
#  DBplus <- NULL
#  DBminus <- NULL
#  stranded <- ifelse(class(reads)=='list', reads[[1]]@stranded, reads@stranded)
#  if(any(strand(DB@islands@unlistData)=='+')) DBplus <- casper:::genomeBystrand(DB, "+")
#  if(any(strand(DB@islands@unlistData)=='-')) DBminus <- casper:::genomeBystrand(DB, "-")
#  if(any(strand(DB@islands@unlistData)=='*')) DBmix <- casper:::genomeBystrand(DB, "*")
#  if(stranded){
#    newexplus <- NULL
#    newexminus <- NULL
#    if(!is.null(DBplus)) newexplus <- findNewExons(casper:::subsetPbam(reads, "+"), DBplus, cov=cov, minConn=minConn, minJunx=minJunx, minLen=minLen, mc.cores=mc.cores)
#    if(!is.null(DBminus)) newexminus <- findNewExons(casper:::subsetPbam(reads, "-"), DBminus, cov=cov, minConn=minConn, minJunx=minJunx, minLen=minLen, mc.cores=mc.cores)
#    if(!is.null(DBmix)) newexmix <- findNewExons(reads, DBmix, cov=cov, minConn=minConn, minJunx=minJunx, minLen=minLen, mc.cores=mc.cores)
#    if(!is.null(newexplus) | !is.null(newexminus)) somex <- 1
#  } else {
#    newex <- findNewExons(reads, DB, cov=cov, minConn=minConn, minJunx=minJunx, minLen=minLen, mc.cores=mc.cores)
#    if(!is.null(newex)) somex <- 1
#  }
#  if(!is.null(somex)){
  #  if(!is.null(DBplus)){
  #    cat("Done...\nCreating denovo genome for positive strand\n")
  #    if(stranded) {
  #      chrs <- as.list(unique(as.character(seqnames(reads@pbam)@values)))
  #      denovoplus <- assignExons2Gene(exons=newexplus, DB=DBplus, reads=casper:::subsetPbam(reads, "+"), chrs=chrs, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
  #    } else {
  #      newexplus <- newex
  #      if(!is.null(DBminus)) {
  #        new <- newex[!(newex %in% DBplus@exonsNI)]
  #        new <- new[!(new %in% DBminus@exonsNI)]
  #        newexplus <- c(newex[newex %in% DBplus@exonsNI], new)
  #      }
  #      chrs <- as.list(unique(as.character(seqnames(reads@pbam)@values)))
  #      denovoplus <- assignExons2Gene(exons=newexplus, DB=DBplus, reads=reads, chrs=chrs, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
  #    }
  #  }
  #  if(!is.null(DBminus)){
  #    cat("Done...\nCreating denovo genome for negative strand\n")
  #    if(stranded) {
  #      chrs <- as.list(unique(as.character(seqnames(reads@pbam)@values)))
  #      denovominus <- assignExons2Gene(exons=newexminus, DB=DBminus, reads=casper:::subsetPbam(reads, "-"), chrs=chrs, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
 #     } else {
 #       newexminus <- newex
 #       if(!is.null(DBplus)) {
 #         new <- newex[!(newex %in% DBminus@exonsNI)]
 #         new <- new[!(new %in% DBplus@exonsNI)]
 #         newexminus <- c(newex[newex %in% DBminus@exonsNI], new)
 #       }
 #       chrs <- as.list(unique(as.character(seqnames(reads@pbam)@values)))
 #       denovominus <- assignExons2Gene(exons=newexminus, DB=DBminus, reads=reads, chrs=chrs, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
 #     }
 #   }
    #if(!is.null(DBminus) & !is.null(DBplus)){
    #  cat("Done...\nMerging denovo genome\n")
    #  browser()
    #  denovo <- mergeStrDenovo(denovoplus, denovominus)
    #} else {
    #  if(!is.null(DBplus)) denovo <- denovoplus
    #  if(!is.null(DBminus)) denovo <- denovominus
    #}
#    chrs <- as.list(unique(as.character(seqnames(reads@pbam)@values)))
#    denovo <- assignExons2Gene(exons=newex, DB=DB, reads=reads, chrs=chrs, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
#  } else {
#    denovo <- DB
#    denovo@denovo <- TRUE
#  }
#  denovo

#}

 
