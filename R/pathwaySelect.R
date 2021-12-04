##' The \code{pathSelect} is function to select pathways of interest for
##' pathway-level congruence analysis by meta pathway enrichement analysis.
##' @title Select pathways of interest for pathway-level congruence analysis
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param pathway.list: list of pathway database. This can be any pathway list from external
##' sources or selected from the contained pathway data in CAMO package.
##' @param pathwaysize.lower.cut: pathway size lower bound cutoff;
##' @param pathwaysize.upper.cut: pathway size upper bound cutoff;
##' @param overlapsize.cut: the lower bound cutoff of overlap size between genes
##' from input data and genes from a pathway.
##' @param med.de.cut: the lower bound cutoff of median number of DE genes in a
##' pathway.
##' @param med.de.cut: the lower bound cutoff of minimum number of DE genes in a
##' pathway.
##' @param qfisher.cut: fisher q-value cutoff from the meta enrichment
##' analysis.
##' @param topPath.indStudy.num: if not NULL, only the union of top pathways is
##' considered.

##' @return a vector of selected pathway names.
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step (see the example in function 'merge')
##' data(human.pathway.list) ## include pathway.list
##' select.pathway = pathSelect(mcmc.merge.list,pathway.list,
##'                             pathwaysize.lower.cut = 5,
##'                             pathwaysize.upper.cut=200,
##'                             overlapsize.cut = 5, med.de.cut =3,
##'                             qfisher.cut = 0.05)
##' select.pathway.list = pathway.list[select.pathway]
##' }


pathSelect <- function(mcmc.merge.list,pathway.list,
                       pathwaysize.lower.cut=10,
                       pathwaysize.upper.cut=200,
                       overlapsize.cut=10,med.de.cut=5,min.de.cut=0,
                       qfisher.cut = 0.05,
                       topPath.indStudy.num = NULL){
  M <- length(mcmc.merge.list)
  summary.list <- vector("list",M)
  for(m in 1:M){
    print(paste("dataset:",m))
    dat <- mcmc.merge.list[[m]]
    deIndex <-  attr(dat,"DEindex")
    summary.list[[m]] <- pathEnrich(dat,deIndex,pathway.list)
  }

  meta.summary <- metaPath(summary.list,pathway.list)

  select.ind <- which(meta.summary$pathway.size >= pathwaysize.lower.cut &
                        meta.summary$pathway.size < pathwaysize.upper.cut &
                        meta.summary$min.overlap.size >=overlapsize.cut &
                        meta.summary$med.DE.inpathway >=med.de.cut &
                        meta.summary$min.DE.inpathway >=min.de.cut &
                        meta.summary$q.fisher < qfisher.cut)

  if(!is.null(topPath.indStudy.num)){
    select.top.pathway <- unique(unlist(lapply(summary.list, function(x){
      subdatai <- x[x$pathway.size<pathwaysize.upper.cut & x$pathway.size>=pathwaysize.lower.cut,]
      sig.pathway <- rownames(subdatai)[order(subdatai$p.value,decreasing=F)[1:topPath.indStudy.num]]
      return(sig.pathway)
    })))
    select.ind = intersect(which(row.names(meta.summary) %in% select.top.pathway),select.ind)
  }

  select.pathways <- rownames(meta.summary)[select.ind]
  return(select.pathways)
}

pathEnrich <- function(dat, deIndex, pathway.list){
  ## dat is data matrix (either rawData or pData) from the individual model for pathway analysis, deIndex is the DE index
  geneDat <- rownames(dat)
  K <- length(pathway.list)
  pathwayName <- names(pathway.list)
  pathwaySize <- sapply(pathway.list,length)

  bk <- intersect(unique(unlist(pathway.list)),geneDat)
  de <- intersect(unique(unlist(pathway.list)),geneDat[deIndex])

  out <- gsa.fisher(x= de, background = bk ,
                    pathway = pathway.list)

  overlapSize <- as.integer(out$DE_in_Set)+as.integer(out$NonDE_in_Set)
  DEnum <- as.integer(out$DE_in_Set)
  OR <- (as.numeric(out$DE_in_Set)*as.numeric(out$NonDE_not_in_Set))/
    (as.numeric(out$DE_not_in_Set)*as.numeric(out$NonDE_in_Set))
  OR <- sapply(OR, function(x) {
    if(is.na(x) || is.infinite(x)) {
      x <- 0 } else{
        x <- x
      }
  }, simplify = T)

  logOR <- ifelse(OR==0,1,log(OR))
  pvalue <- as.numeric(out$pvalue)

  summary <- data.frame(pathway.size= pathwaySize, overlap.size= overlapSize,
                        DE.inpathway = DEnum, Odds.ratio = OR,
                        log.odds.ratio = logOR, p.value=pvalue)
  rownames(summary) <- pathwayName
  return(summary)
}


metaPath <- function(pathSummary, pathway.list){

  M <- length(pathSummary)
  K <- length(pathway.list)
  pathwayName <- names(pathway.list)
  pathwaySize <- sapply(pathway.list,length)

  pdat <- Reduce("cbind",lapply(pathSummary,function(x) x$p.value))
  rownames(pdat) <- pathwayName
  p.fisher <- apply(pdat,1,fisher)
  q.fisher <- p.adjust(p.fisher,method="BH")

  overlapMinSize <- apply(Reduce("cbind",lapply(pathSummary,
                                                function(x) x$overlap.size)),
                          1,min)

  overlapMedDE <-  apply(Reduce("cbind",lapply(pathSummary,
                                               function(x) x$DE.inpathway)),1,
                         function(x) (
                           if(M %% 2 ==0){
                             (sort(x)[M/2]+sort(x)[(M/2+1)])/2
                           } else { sort(x)[(M+1)/2] } ))

  overlapMinDE <- apply(Reduce("cbind", lapply(pathSummary,
                                               function(x) x$DE.inpathway)), 1, min)

  meta.summary <- data.frame(pathway.size= pathwaySize,
                             min.overlap.size= overlapMinSize,
                             med.DE.inpathway = overlapMedDE,
                             min.DE.inpathway = overlapMinDE,
                             p.fisher=as.numeric(p.fisher),
                             q.fisher=as.numeric(q.fisher))
  rownames(meta.summary) <- pathwayName
  return(meta.summary)

}



