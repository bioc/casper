simReads <- function(islandid, nSimReads, pis, rl, seed, writeBam, distrs, genomeDB, repSims=FALSE, bamFile=NULL, stranded=FALSE, verbose=TRUE, chr=NULL, mc.cores=1){
  ##### Simulate reads
  if(writeBam==1) {
    chrlen <- seqlengths(genomeDB@exonsNI)
    tmp <- sub(" ", ".", Sys.time())
    if(is.null(chr)) lr_file <- paste(bamFile, tmp, ".sam", sep="")
  }
  if(!is.null(chr)){
    isl2chr <- genomeDB@exon2island$seqname
    names(isl2chr) <- genomeDB@exon2island$island
    lr_file=NULL
    sims <- lapply(chr, function(x){
      if(verbose) cat("Simulating ", x, "\n")
      if(writeBam) lr_file <- paste(bamFile, tmp, ".", x, ".sam", sep="")
      islandid <- islandid[islandid %in% names(isl2chr)[isl2chr %in% x]]
      if(writeBam) sims <- casperSim(genomeDB=genomeDB, distrs=distrs, nSimReads=nSimReads, pis=pis, islandid=islandid, lr_file=lr_file, rl=rl, chrlen=chrlen, seed=seed, bam=writeBam, chr=x, verbose=verbose)
    if(!writeBam) sims <- casperSim(genomeDB=genomeDB, distrs=distrs, nSimReads=nSimReads, pis=pis, islandid=islandid, rl=rl, seed=seed, bam=writeBam, chr=x, verbose=verbose)
    })
    sims <- do.call('c', sims)
  } else {
    if(writeBam) {
      sims <- casperSim(genomeDB=genomeDB, distrs=distrs, nSimReads=nSimReads, pis=pis, islandid=islandid, lr_file=lr_file, rl=rl, chrlen, seed, bam=writeBam, verbose=verbose)
    } else sims <- casperSim(genomeDB=genomeDB, distrs=distrs, nSimReads=nSimReads, pis=pis, islandid=islandid, rl=rl, seed=seed, bam=writeBam,  verbose=verbose)
  }
    if (verbose) cat("Splitting counts\n")
    counts <- splitPaths(sims$pc, genomeDB, mc.cores=mc.cores, stranded=stranded, geneid=islandid)
    pc <- new("pathCounts", counts=counts, denovo=FALSE, stranded=stranded)
   
  if(repSims==1) {
    sims$pc <- NULL
    sims$Nsim <- NULL
    ans <- list(sims=sims, pc=pc)
  }
  else ans <- pc
  ans
}


splitPaths <- function(paths, DB, mc.cores, stranded, geneid){
    paths <- paths[grepl("^\\..*\\.$", names(paths))]
    sel <- strsplit(names(paths), split='-|\\.')
    sel <- lapply(sel, function(x) x[x!=''])
    sel1 <- lapply(sel, "[", 2)
    sel1 <- unlist(sel1)
    islands <- DB@islands[geneid]
    nislEx <- elementLengths(islands)
    nislEx <- rep(names(islands), nislEx)
    islEx <- names(islands@unlistData)
    names(islEx) <- nislEx
    isl <- match(sel1, islEx)
    isl <- names(islEx)[isl]
    splCounts <- split(paths, isl)
    splCounts <- lapply(splCounts, function(x) x[grepl("-", names(x))])
  if(DB@denovo){
    sel <- lapply(sel, "[", -1)
    tmp <- split(sel, isl)
    if(mc.cores>1) {
      require(parallel)
      tmp1 <- parallel::mclapply(names(tmp), function(x){
        n <- sapply(tmp[[x]], length)
        nn <- unlist(tmp[[x]])
        names(nn) <- rep(names(splCounts[[x]]), n)
        nnn <- nn %in% names(DB@islands[[x]])
        nnnn <- tapply(nnn, names(nn), all)
        splCounts[[x]][names(splCounts[[x]]) %in% names(nnnn)[nnnn]]
      }, mc.cores=mc.cores)
    } else {
      tmp1 <- lapply(names(tmp), function(x){
        n <- sapply(tmp[[x]], length)
        nn <- unlist(tmp[[x]])
        names(nn) <- rep(names(splCounts[[x]]), n)
        nnn <- nn %in% names(DB@islands[[x]])
        nnnn <- tapply(nnn, names(nn), all)
        splCounts[[x]][names(splCounts[[x]]) %in% names(nnnn)[nnnn]]
      })
    }
    names(tmp1) <- names(tmp)
    splCounts <- tmp1
  }
  ans <- vector(length(islands), mode='list')
  names(ans) <- names(islands)
  ans[names(splCounts)] <- splCounts
  if(!stranded) ans <- list(ans)
  else {
    is <- as.character(strand(DB@islands@unlistData))[cumsum(c(1, elementLengths(DB@islands)[-length(DB@islands)]))]
    names(is) <- names(DB@islands)    
    plus <- ans[names(ans) %in% names(DB@islands)[is=='+']]
    minus <- ans[names(ans) %in% names(DB@islandStrand)[is=='-']]
    ans <- list(minus=minus, plus=plus)
  }
  ans
}

casperSim <- function(genomeDB, distrs, nSimReads, pis, islandid, lr_file=NULL, rl, chrlen, seed, bam, chr=NULL, verbose=FALSE){
  if (verbose) cat("Formatting input\n")

  sel <- islandid
  nSimReads <- nSimReads[sel]
  txs <- genomeDB@transcripts[sel]
  sel1 <- unlist(lapply(txs, is.null))
  txs <- txs[!sel1]
  variant_num <- unlist(lapply(txs, length))
  starts <- start(genomeDB@exonsNI)
  names(starts) <- names(genomeDB@exonsNI)
  ends <- end(genomeDB@exonsNI)
  names(ends) <- names(starts)
  chroms <- as.character(seqnames(genomeDB@exonsNI))
  names(chroms) <- names(genomeDB@exonsNI)
  exon_num <- unlist(lapply(txs, function(x) lapply(x, length)))
  exs <- unlist(txs)
  exon_st <- starts[as.character(exs)]
  exon_end <- ends[as.character(exs)]
  exs <- unname(exs)
  sel2 <- unlist(lapply(lapply(txs, '[[', 1), '[', 1))
  chroms <- as.character(chroms[as.character(sel2)])
  tmps <- rep(1:length(exon_num), exon_num)
  widths <- exon_end-exon_st+1
  variant_len <- tapply(widths, tmps, sum)
  tmptxs <- unlist(txs, recursive=F)
  txs <- unlist(lapply(txs,names))
  tx_strand <- as.integer(sapply(tmptxs, function(x) ifelse(x[1]<x[2], 1, -1)))
  if(!all(txs %in% names(pis))) stop("Wrong pis vector, some transcripts missing")
  if(!all(nSimReads>0)) stop("nSimReads with zero entries")
  ge <- nSimReads[sel]
  ge <- as.integer(rep(0:(length(ge)-1), ge))
  ve=pis[txs]
  vn=variant_num
  vl=as.integer(variant_len)
  en=exon_num
  es=exon_st
  ee=exon_end
  ei=exs
  ngenes=length(sel)
  #ldv <- sample(as.integer(names(distrs@lenDis)), prob=distrs@lenDis/sum(distrs@lenDis), size=sum(nSimReads), replace=T)
  #ldd <- as.numeric(1
  #browser()
  ldv <- cumsum(distrs@lenDis)/sum(distrs@lenDis)
  #th <- seq(0,1, length=length(ldv))
  #ldv <- approxfun(ldv, ldd)
  #th <- seq(0, 1, len=10000)
  #ldv <- ldv(th)
  #ldv[1] <- 0
  ldd <- as.numeric(names(distrs@lenDis))
  #ldd <- round(seq(min(ldd), max(ldd), length=10000))
  th <- seq(0, 1, len=10000)
  std <- distrs@stDis(th)
  std[1] <- 0
  std[length(th)] <- 1
  sdv <- approxfun(std[!is.na(std)], th[!is.na(std)])(th)
  sdd <- std
  if(bam) write.sam.header(chrlen[names(chrlen) %in% unique(chroms)], lr_file, max(vn))
  if (verbose) cat("Simulating fragments\n")
# insideBam deprecated, only in case it is useful in simulations, return from C all information written to the bam file
  insideBam=as.integer(0)
  #browser()

  #dyn.load("/Volumes/biostats/cstephan/casper_bioC/casper/src/casper.so")
  #source("/Volumes/biostats/cstephan/casper_bioC/casper/R/simPost.R")
  ans <- .Call("casperSimC", ge, ve, vn, as.integer(vl), en, es, ee, ei, ldv, ldd, sdv, sdd, as.integer(rl), length(ge), as.integer(tx_strand), lr_file, chroms, as.integer(seed), as.integer(bam), as.integer(insideBam), as.integer(verbose))

#pdf("try.dis.pdf");plot(density(ans[[3]]), xlim=c(0,500));lines(density(ans2[[3]]), col=2);lines(density(ans3[[3]]), col=3);dev.off()
  
  ans <- ans[1:7]
  names(ans[[7]]) <- ans[[6]]
  ans[[6]] <- NULL
  names(ans) <- c("varl", "st", "len", "abst", "strand", "pc")
  ans
}

write.sam.header <- function(chrlen, file, nvars){
  chr.line <- function(chr, len, file, append) cat("@SQ\tSN:", chr, "\tLN:", len, "\n", sep="", file=file, append=append)
  chr.line(names(chrlen)[1], chrlen[1], file, append=F)
  if(length(chrlen)>1) for(i in 2:length(chrlen)) chr.line(names(chrlen)[i], chrlen[i], file, append=T)
  cat("@PG\tID:simulation\tPN:simulation\tVN:0.0.1\n", file=file, append=T)
  for(i in 1:nvars) cat("@RG\tID:", i, "\tSM:", i, "\n", file=file, append=T)
}

