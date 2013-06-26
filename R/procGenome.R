require(methods)

## METHODS

setMethod("getIsland",signature(entrezid='character',txid='missing',genomeDB='annotatedGenome'),function(entrezid, txid, genomeDB) {
  as.character(genomeDB@aliases[genomeDB@aliases$gene_id==entrezid,'island_id'][1])
}
)

setMethod("getIsland",signature(entrezid='missing',txid='character',genomeDB='annotatedGenome'),function(entrezid, txid, genomeDB) {
  as.character(genomeDB@aliases[txid,'island_id'])
}
)

setMethod("getChr",signature(entrezid='missing',txid='missing',islandid='character', genomeDB='annotatedGenome'), function(entrezid, txid, islandid, genomeDB) {
  as.character(genomeDB@exon2island$seqnames[match(TRUE,genomeDB@exon2island$island == islandid)])
}
)

setMethod("getChr",signature(entrezid='character',txid='missing',islandid='missing', genomeDB='annotatedGenome'), function(entrezid, txid, islandid, genomeDB) {
  islandid <- getIsland(entrezid=entrezid, genomeDB=genomeDB)
  getChr(islandid=islandid, genomeDB=genomeDB)
}
)

setMethod("getChr",signature(entrezid='missing',txid='character',islandid='missing', genomeDB='annotatedGenome'), function(entrezid, txid, islandid, genomeDB) {
  islandid <- getIsland(txid=txid, genomeDB=genomeDB)
  getChr(islandid=islandid, genomeDB=genomeDB)
}
)

setMethod("transcripts", signature(entrezid='missing',islandid='character',genomeDB='annotatedGenome'), function(entrezid, islandid, genomeDB) {
  IRangesList(lapply(genomeDB@transcripts[[islandid]],function(z) ranges(genomeDB@islands[[islandid]][as.character(z)])))
}
)

setMethod("transcripts", signature(entrezid='character',islandid='missing',genomeDB='annotatedGenome'), function(entrezid, islandid, genomeDB) {
  txids <- rownames(genomeDB@aliases)[genomeDB@aliases$gene_id==entrezid]
  islandid <- getIsland(entrezid=entrezid,genomeDB=genomeDB)
  IRangesList(lapply(genomeDB@transcripts[[islandid]],function(z) ranges(genomeDB@islands[[islandid]][as.character(z)])))[txids]
}
)

setGeneric("txLength", function(islandid, txid, genomeDB) standardGeneric("txLength"))
setMethod("txLength", signature(islandid='missing', txid='missing', genomeDB='annotatedGenome'), function(islandid, txid, genomeDB) {
  txid <- rownames(genomeDB@aliases)
  txLength(txid=txid, genomeDB=genomeDB)
}
)

setMethod("txLength", signature(islandid='character', txid='missing', genomeDB='annotatedGenome'), function(islandid, txid, genomeDB) {
  txid <- rownames(genomeDB@aliases)[genomeDB@aliases$island_id==islandid]
  txLength(txid=txid, genomeDB=genomeDB)
}
)

setMethod("txLength", signature(islandid='missing', txid='character', genomeDB='annotatedGenome'), function(islandid, txid, genomeDB) {
  if(length(genomeDB@txLength)==0){
    tx <- unlist(genomeDB@transcripts,recursive=FALSE)
    names(tx) <- sapply(strsplit(names(tx),'\\.'),'[[',2)
    tx <- tx[txid]
    tx <- data.frame(tx=rep(names(tx),sapply(tx,length)), exon=unlist(tx))
    #Get exon id & width into data.frame
    w <- unlist(width(genomeDB@islands))
    r <- unlist(ranges(genomeDB@islands))
    names(w) <- sapply(strsplit(names(r),'\\.'),'[[',2)
    w <- data.frame(exon=as.integer(names(w)), width=w)
    #Merge and by
    exonw <- merge(tx,w,by='exon',all.x=TRUE)
    txw <- by(exonw$width, INDICES=factor(exonw$tx), FUN=sum)
    ans <- as.integer(txw); names(ans) <- names(txw)
  } else ans <- genomeDB@txLength[txid]
  return(ans)
}
)

setGeneric("subsetGenome", function(islands, chr, genomeDB) standardGeneric("subsetGenome"))
setMethod("subsetGenome", signature(islands='character', chr='missing', genomeDB='annotatedGenome'), function(islands, chr, genomeDB) {
  islands <- unique(islands)
  isl <- subset(genomeDB@islands, names(genomeDB@islands) %in% islands)
  txs <- unlist(sapply(genomeDB@transcripts[islands], names))
  exs <- unique(names(genomeDB@islands[islands]@unlistData))
  alia <- genomeDB@aliases[txs,]
  ex2is <- genomeDB@exon2island[as.character(genomeDB@exon2island$island) %in% islands,]
  txs <- genomeDB@transcripts[islands]
  exs <- subset(genomeDB@exonsNI, names(genomeDB@exonsNI) %in% exs)
  new("annotatedGenome", islands=isl, transcripts=txs, exonsNI=exs, aliases=alia, exon2island=ex2is, dateCreated=genomeDB@dateCreated, denovo=genomeDB@denovo, genomeVersion=genomeDB@genomeVersion)
})
setMethod("subsetGenome", signature(islands='missing', chr='character', genomeDB='annotatedGenome'), function(islands, chr, genomeDB) {
  islands <- unique(as.character(genomeDB@exon2island[genomeDB@exon2island$seqnames %in% chr,]$island))
  subsetGenome(islands=islands, genomeDB=genomeDB)
})

### FUNCTIONS

makeIslands <- function(exons){
  txs <- as.integer(as.factor(names(exons)))
  totEx <- length(exons)
  uniex <- unique(exons)
  nexR <- length(uniex)
  islands <- rep(0, nexR)
  nexons <- names(exons)
  exons <- as.character(sprintf("%d",exons))
  names(exons) <- nexons
  tabex <- table(exons)
  tabex <- tabex[exons]
  tabtx <- table(names(exons))
  tabtx <- tabtx[names(exons)]
  ans<-.Call("makeGeneIslands", as.integer(exons), islands, uniex, txs, totEx, nexR, as.integer(tabex), as.integer(tabtx))
  names(ans) <- as.character(sprintf("%d",uniex))
  ans
}

generateNOexons<-function(exByTx, startId=1, mc.cores){
  exByTx<-unlist(exByTx)
  strand(exByTx) <- "*"
  exonsNI <- disjoin(exByTx)
  names(exonsNI) <- startId:(length(exonsNI)+startId-1)
  overEx <- findOverlaps(exonsNI, exByTx)
  exid <- names(exonsNI)[queryHits(overEx)]
  exkey <- split(exid, values(exByTx)$exon_id[subjectHits(overEx)])
  if(mc.cores>1) require(parallel)
  if(mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) {
      exkey<-parallel::mclapply(exkey, unique, mc.cores=mc.cores)
    } else stop('parallel library has not been loaded!')
  }
  else {
    exkey<-lapply(exkey, unique, mc.cores=mc.cores)
  }
  res<-list(exkey=exkey, exons=exonsNI)
  res 
}

genomeBystrand <- function(DB, strand){
  is <- as.character(strand(DB@islands@unlistData))[cumsum(c(1, elementLengths(DB@islands)[-length(DB@islands)]))]
  sel <- names(DB@islands)[is==strand]
  islands <- DB@islands[sel]
  transcripts <- DB@transcripts[sel]
  exonsNI <- DB@exonsNI[names(DB@exonsNI) %in% as.character(unlist(transcripts)),]
  exon2island <- DB@exon2island[rownames(DB@exon2island) %in% as.character(unlist(transcripts)),]
  txid <- unlist(lapply(transcripts, names))
  aliases <- DB@aliases[DB@aliases$tx %in% txid,]
  ans <- new("annotatedGenome", aliases=aliases, denovo=TRUE, exonsNI=exonsNI, transcripts=transcripts, exon2island=exon2island, dateCreated=Sys.Date(), genomeVersion=DB@genomeVersion, islands=islands)
  ans
}


rmDuplicateTxs <- function(txs, txname, exid) {
  txs <- txs[match(unique(mcols(txs)[,txname]), mcols(txs)[,txname]),]
  exid <- exid[mcols(txs)[,txname]]
  txs <- txs[match(unique(exid), exid),]
  aliases <- values(txs)[1:3]
  aliases[,3] <- unlist(aliases[,3])
  aliases <- as.data.frame(aliases, stringsAsFactors=FALSE)
  aliases <- cbind(aliases, exid=exid[as.character(aliases[,txname])])
  alitx <- split(as.character(aliases[,2]), aliases[,4])
  names(alitx) <- unlist(lapply(alitx, function(x) ifelse(length(x)>1, x[x %in% txs@elementMetadata[,txname]], x)))
  nalitx <- rep(names(alitx), unlist(lapply(alitx, length)))
  alitx <- unlist(alitx)
  names(alitx) <- nalitx
  aliases <- cbind(aliases, tx=names(alitx)[match(aliases[,2], alitx)])
  return(list(txs=txs,aliases=aliases))
}


createGenome <- function(txs, Exons, genome, mc.cores) {
  exid <- sapply(txs@elementMetadata$exon_id, function(x) paste(unlist(x), collapse="."))
  names(exid) <- values(txs)[,'tx_name']
  txs <- rmDuplicateTxs(txs, txname='tx_name', exid=exid)
  aliases <- txs$aliases; txs <- txs$txs
  Exons<-Exons[names(Exons) %in% txs@elementMetadata$tx_id,]
  txnames<-match(names(Exons), txs@elementMetadata$tx_id)
  names(Exons)<-txs@elementMetadata$tx_name[txnames]
  cat("Finding non-overlapping exons\n")
  txStrand <- as.character(strand(txs))
  names(txStrand) <- mcols(txs)[,"tx_name"] #values(txs)$tx_name
  exonsNI <- generateNOexons(Exons, startId=1, mc.cores=mc.cores)  
  exkey <- exonsNI$exkey
  exonsNI <- exonsNI$exons
  #Find transcript structure for new exons
  cat("Remapping transcript structure to new exons\n")
  s<-ifelse(as.character(Exons@unlistData@strand)=="+", 1, -1)
  exids <- values(Exons@unlistData)$exon_id
  idx <- values(Exons@unlistData)$exon_rank
  idx <- cumsum(idx==1)
  exids <- s*exids
  exon_ids <- tapply(exids, idx, function(x) if(any(x<0)) {rev(-x)} else x)
  names(exon_ids) <- names(Exons)
  exlen <- sapply(exkey, length)
  newTxs <- as.integer(unlist(exkey[as.character(sprintf("%d",unlist(exon_ids)))]))
  lenkey <- exlen[as.character(sprintf("%d",unlist(exon_ids)))]
  lenkey <- tapply(lenkey, rep(names(exon_ids), sapply(exon_ids, length)), sum)
  lenkey <- lenkey[names(exon_ids)]
  newTxs <- split(newTxs, rep(names(lenkey), lenkey))
  newTxs <- newTxs[names(Exons)]
  geneids <- mcols(txs)[,"gene_id"] #values(txs)$gene_id
  #geneids <- unlist(geneids)
  names(geneids)<- mcols(txs)[,"tx_id"] #unlist(values(txs)$tx_id)
  # Make islands
  ex2tx <- unlist(newTxs)
  names(ex2tx) <- rep(names(newTxs), unlist(lapply(newTxs, length)))
  islands <- casper:::makeIslands(ex2tx)
  cat("Splitting transcripts\n")
  extxs <- unlist(lapply(newTxs, "[", 1))    
  sel <- match(extxs, names(exonsNI))
  tx2island <- islands[as.character(sel)]
  names(tx2island) <- names(newTxs)
  transcripts <- newTxs[names(tx2island)]
  sel <- txStrand[names(transcripts)]=='-'
  transcripts[sel] <- sapply(transcripts[sel], rev)
  transcripts <- split(transcripts, tx2island)
  islandStrand <- tapply(txStrand, tx2island[names(txStrand)], unique)
  ss <- sapply(islandStrand, length)
  islandStrand[ss>1] <- "*"
  id2tx <- data.frame(island_id=rep(names(transcripts),sapply(transcripts,length)) , txname=unlist(sapply(transcripts,names)))
  rownames(id2tx) <- id2tx$txname
  aliases$island_id <- id2tx[rownames(aliases),'island_id']
  exon2island <- as.data.frame(ranges(exonsNI))
  exon2island <- exon2island[,-4]
  exon2island$seqnames <- as.character(seqnames(exonsNI))
  rownames(exon2island) <- names(exonsNI)
  exon2island$island <- islands[rownames(exon2island)]
  exp <- exonsNI[unlist(islandStrand[as.character(exon2island$island)]) == '+']
  strand(exp) <- "+"
  tmp <- split(exp, islands[names(exp)])
  exu <- exonsNI[unlist(islandStrand[as.character(exon2island$island)]) == '*']
  strand(exu) <- "*"
  tmpu <- split(exu, islands[names(exu)])
  exm <- exonsNI[islandStrand[as.character(exon2island$island)]=='-']
  strand(exm) <- "-"
  tmpm <- split(rev(exm), islands[names(rev(exm))])
  tmp <- c(tmp, tmpm, tmpu)
  tmp <- tmp[names(transcripts)]
  ans <- new("annotatedGenome", islands=tmp, transcripts=transcripts, exon2island=exon2island, aliases=aliases, exonsNI=exonsNI, dateCreated=Sys.Date(), genomeVersion=genome, denovo=FALSE)
  txL <- txLength(genomeDB=ans)
  ans@txLength <- txL
  return(ans)
}


setMethod("procGenome", signature(genDB="TranscriptDb"), function(genDB, genome, mc.cores=1) {
#  genDB<-makeTranscriptDbFromUCSC(genome=genome, tablename="refGene")
  cat("Processing Exons and Transcrips\n")
  txs <- GenomicFeatures::transcripts(genDB,columns=c("tx_id","tx_name","gene_id","exon_id"))
  Exons <- exonsBy(genDB, by="tx")
  createGenome(txs=txs, Exons=Exons, genome=genome, mc.cores=mc.cores)
 } 
)


setMethod("procGenome", signature(genDB="GRanges"), function(genDB, genome, mc.cores=1) {
  cat("Formatting GTF table with GenomicFeatures tools...\n")
  if (all(c("gene_id", "transcript_id") %in% colnames(mcols(genDB)))) {
    tables <- GenomicFeatures:::.prepareGTFTables(genDB, exonRankAttributeName=NULL)
  } else {
    stop("Columns named 'gene_id' and 'transcript_id' not found")
  }
  chroms <- unique(tables[["transcripts"]][["tx_chrom"]])
  chrominfo <- data.frame(chrom = chroms, length = rep(NA,length(chroms)))
  #Eliminate transcripts with unknown strands (usually transcripts with 1 exon)
  sel <- (tables$transcripts$tx_strand %in% c('+','-'))
  tables$transcripts <- tables$transcripts[sel,]
  sel <- (tables$splicings$exon_strand %in% c('+','-'))
  tables$splicings <- tables$splicings[sel,]
  tables$genes <- tables$genes[tables$genes$tx_name %in% tables$transcripts$tx_name,]
  cat("Making TranscriptDb object...\n")
  txdb <- makeTranscriptDb(transcripts=tables[["transcripts"]], splicings=tables[["splicings"]], genes=tables[["genes"]], chrominfo=chrominfo, reassign.ids=TRUE)
  txs <- GenomicFeatures::transcripts(txdb,columns=c("tx_id","tx_name","gene_id","exon_id"))
  Exons <- exonsBy(txdb, by="tx")
  ans <- createGenome(txs=txs, Exons=Exons, genome=genome, mc.cores=mc.cores)
  return(ans)
}
)

