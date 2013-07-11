splitGenomeByLength <- function(DB, breaks=c(0,3000,5000,Inf)) {
  islands <- as.data.frame(DB@islands)
  islandw <- sqldf::sqldf("select element, sum(width) from islands group by element")
  #if (missing(breaks)) breaks <- c(0,quantile(islandw[,2],probs=c(1/3,2/3)),Inf)
  l <- cut(islandw[,2], breaks=breaks)
  names(l) <- islandw[,1]
  l <- l[names(DB@islands)]
  ans <- vector("list",length(levels(l))); names(ans) <- as.character(levels(l))
  for (i in 1:length(ans)) {
    sel <- names(l)[l==levels(l)[i]]
    newDB <- DB
    newDB@islands <- newDB@islands[sel]
    newDB@transcripts <- newDB@transcripts[sel]
    newDB@exon2island <- newDB@exon2island[newDB@exon2island$island %in% as.integer(sel),]
    newDB@exonsNI <- newDB@exonsNI[newDB@exonsNI %over% newDB@islands,]
    newDB@aliases <- newDB@aliases[as.integer(newDB@aliases$island_id) %in% as.integer(sel),]
    newDB@txLength <- newDB@txLength[names(newDB@txLength) %in% unlist(sapply(newDB@transcripts, names))]
    ans[[i]] <- newDB
  }
  return(ans)
}
