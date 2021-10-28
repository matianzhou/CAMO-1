##' The \code{multi_ACS_ADS_pathway} is function to perform congruence analysis
##' for multiple pairs, generating pathway-specific c-scores, d-scores and their permuted p-value.
##' @title Pathway level congruence analysis for multiple pairs
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param dataset.names: a vector of dataset names.
##' @param measure: three types of scores to be used: "youden",
##' "Fmeasure","geo.mean". Default is "Fmeasure".
##' @param B: number of permutations.
##' @param parallel: whether to perform parallel ('mclapply' will be used if parallel=T) computing in permutation.
##' @param n.cores: if parallel=T, the number of cores to use (mc.cores parameter in 'mclapply' function).
##' @return Four lists:
##' \itemize{
##' \item ACS.mat: pathway-specific c-scores.
##' \item ACSpvalue.mat: p-values for pathway-specific c-scores.
##' \item ADS.mat: pathway-specific d-scores.
##' \item ADSpvalue.mat: p-values for pathway-specific d-scores.
##' }
##' In addition, the four data matrices are written to a folder named "ACS_ADS_Pathway".
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step (see the example in function 'merge')
##' #select.pathway.list from the pathSelect step (see the example in function 'pathSelect')
##' dataset.names = c("hb","hs","ht","ha","hi","hl",
##'                   "mb","ms","mt","ma","mi","ml")
##' ACS_ADS_pathway = multi_ACS_ADS_pathway(mcmc.merge.list,dataset.names,
##'                                         select.pathway.list,
##'                                         measure="Fmeasure",
##'                                         B=1000,parallel=F)
##' }
multi_ACS_ADS_pathway <- function(mcmc.merge.list,dataset.names,
                                  select.pathway.list,
                                  measure="Fmeasure",
                                  B=100,parallel=F,n.cores=4){

  names(mcmc.merge.list) <- dataset.names
  M <- length(mcmc.merge.list)
  P <- choose(M,2)
  select.pathways <- names(select.pathway.list)
  data_genes <- rownames(mcmc.merge.list[[1]])
  pathway.size <- sapply(select.pathway.list,function(x) {
    length(intersect(data_genes,x))})
  K <- length(select.pathways)

  ACS <- ACSpvalue <- ADS <- ADSpvalue <- array(1,dim=c(K,M,M),dimnames=
                                                  list(select.pathways,dataset.names,dataset.names))
  for(i in 1:(M-1)){
    for(j in (i+1):M){
      print(paste("pair: dataset ",i," and dataset ",j,sep=""))

      dat1 <- mcmc.merge.list[[i]]
      dat2 <- mcmc.merge.list[[j]]
      deIndex1 <- attr(dat1,"DEindex")
      deIndex2 <- attr(dat2,"DEindex")
      names(deIndex1) <- rownames(dat1)[deIndex1]
      names(deIndex2) <- rownames(dat2)[deIndex2]
      data_genes <- rownames(dat1)

      print("Calculate ECS & EDS for each pathway...")
      marginOut <- margin_pathway(dat1,dat2,select.pathway.list,
                                  measure=measure)
      print("Permutate genes for each pathway...")
      permOut_pathway <- perm_pathway(dat1,dat2,select.pathway.list,
                                      measure=measure,B=B,parallel=parallel,n.cores=n.cores)
      ACS[,i,j] <- ACS[,j,i] <- acs_pathway <- ACS_pathway(dat1,dat2,deIndex1,deIndex2,
                                                           select.pathway.list,
                                                           measure=measure,marginOut=marginOut)

      ACSpvalue[,i,j] <- ACSpvalue[,j,i] <- pACS_pathway(dat1,dat2,deIndex1,deIndex2,
                                                         select.pathway.list,
                                                         measure=measure,acs=acs_pathway,
                                                         permOut=permOut_pathway,
                                                         marginOut=marginOut)

      ADS[,i,j] <- ADS[,j,i] <- ads_pathway <- ADS_pathway(dat1,dat2,deIndex1,deIndex2,
                                                           select.pathway.list,
                                                           measure=measure,marginOut=marginOut)

      ADSpvalue[,i,j] <- ADSpvalue[,j,i] <- pADS_pathway(dat1,dat2,deIndex1,deIndex2,
                                                         select.pathway.list,
                                                         measure=measure,ads=ads_pathway,
                                                         permOut=permOut_pathway,
                                                         marginOut=marginOut)
    }
  }

  ACS.mat <- ACSpvalue.mat <- ADS.mat <- ADSpvalue.mat <- data.frame(matrix(0,K,P))
  rownames(ACS.mat) <- rownames(ACSpvalue.mat) <- rownames(ADS.mat) <- rownames(ADSpvalue.mat) <- select.pathways

  combinations <- combn(dataset.names,m=2)
  colnames(ACS.mat) <- colnames(ACSpvalue.mat) <- colnames(ADS.mat) <- colnames(ADSpvalue.mat) <- apply(combinations,2,FUN=function(x) paste(x,collapse ="_"))


  for(p in 1:P){
    pairname <- colnames(ACS.mat)[p]
    name1 <- strsplit(pairname,split="_",fixed=T)[[1]][1]
    name2 <- strsplit(pairname,split="_",fixed=T)[[1]][2]
    ACS.mat[,pairname] <- ACS[,name1,name2]
    ACSpvalue.mat[,pairname] <- ACSpvalue[,name1,name2]
    ADS.mat[,pairname] <- ADS[,name1,name2]
    ADSpvalue.mat[,pairname] <- ADSpvalue[,name1,name2]

  }

  dir.path <- "ACS_ADS_Pathway"
  if (!file.exists(dir.path)) dir.create(dir.path)
  write.csv(ACS.mat,file=paste(paste(dir.path,"ACS_pathway_",sep="/"),M,".csv",sep=""))
  write.csv(ACSpvalue.mat,file=paste(paste(dir.path,"ACSpvalue_pathway_",sep="/"),M,".csv",sep=""))
  write.csv(ADS.mat,file=paste(paste(dir.path,"ADS_pathway_",sep="/"),M,".csv",sep=""))
  write.csv(ADSpvalue.mat,file=paste(paste(dir.path,"ADSpvalue_pathway_",sep="/"),M,".csv",sep=""))

  out <- list(ACS.mat=ACS.mat,ACSpvalue.mat=ACSpvalue.mat,ADS.mat=ADS.mat,ADSpvalue.mat=ADSpvalue.mat)
  return(out)
}
