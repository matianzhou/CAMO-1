##' The \code{multi_ACS_ADS_global} is function to perform congruence analysis
##' for multiple study pairs, generating genome-wide c-scores, d-scores and their permuted p-values.
##' @title Genome-wide congruence analysis for multiple pairs
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param dataset.names: a vector of dataset names.
##' @param measure: three types of scores to be used: "youden",
##' "Fmeasure","geo.mean". Default is "Fmeasure".
##' @param B: number of permutations.
##' @return Four lists:
##' \itemize{
##' \item ACS: genome-wide c-scores.
##' \item ACSpvalue: p-values for genome-wide c-scores.
##' \item ADS: genome-wide d-scores.
##' \item ADSpvalue: p-values for genome-wide d-scores.
##' }
##' In addition, the four data matrices are written to a folder named "ACS_ADS_Global".
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step (see the example in function 'merge')
##' dataset.names = c("hb","hs","ht","ha","hi","hl",
##'                   "mb","ms","mt","ma","mi","ml")
##' ACS_ADS_global = multi_ACS_ADS_global(mcmc.merge.list,dataset.names,
##' measure="Fmeasure",B=1000)
##' }

multi_ACS_ADS_global <- function(mcmc.merge.list,dataset.names,
                                 measure="Fmeasure",B=100){

  names(mcmc.merge.list) <- dataset.names
  M <- length(mcmc.merge.list)
  P <- choose(M,2)
  ACS <- ACSpvalue <- matrix(NA,M,M)
  rownames(ACS) <- colnames(ACS) <-
    rownames(ACSpvalue) <- colnames(ACSpvalue) <- dataset.names
  ADS <- ADSpvalue <- matrix(NA,M,M)
  rownames(ADS) <- colnames(ADS) <-
    rownames(ADSpvalue) <- colnames(ADSpvalue) <- dataset.names

  diag(ACS) <- diag(ACSpvalue) <- diag(ADS) <- diag(ADSpvalue)  <- 1

  for(i in 1:(M-1)){
    for(j in (i+1):M){
      dat1 <- mcmc.merge.list[[i]]
      dat2 <- mcmc.merge.list[[j]]
      deIndex1 <- attr(dat1,"DEindex")
      deIndex2 <- attr(dat2,"DEindex")

      permOut <- perm_global(dat1,dat2,measure="Fmeasure",B=B)
      ACS[i,j] <- ACS[j,i] <- acs_global <- ACS_global(dat1,dat2,deIndex1,deIndex2,
                                                       measure=measure)
      ACSpvalue[i,j] <- ACSpvalue[j,i] <- pACS_global(dat1,dat2,deIndex1,deIndex2,
                                                      measure=measure,acs=acs_global,permOut=permOut)
      ADS[i,j] <- ADS[j,i] <- ads_global <- ADS_global(dat1,dat2,deIndex1,deIndex2,
                                                       measure=measure)
      ADSpvalue[i,j] <- ADSpvalue[j,i] <- pADS_global(dat1,dat2,deIndex1,deIndex2,
                                                      measure=measure,ads=ads_global,permOut=permOut)
      print(paste("pair: dataset ",i," and dataset ",j,sep=""))
    }
  }
  dir.path <- "ACS_ADS_Global"
  if (!file.exists(dir.path)) dir.create(dir.path)
  write.csv(ACS,file=paste(paste(dir.path,"ACS_global_",sep="/"),M,".csv",sep=""))
  write.csv(ACSpvalue,file=paste(paste(dir.path,"ACSpvalue_global_",sep="/"),M,".csv",sep=""))
  write.csv(ADS,file=paste(paste(dir.path,"ADS_global_",sep="/"),M,".csv",sep=""))
  write.csv(ADSpvalue,file=paste(paste(dir.path,"ADSpvalue_global_",sep="/"),M,".csv",sep=""))

  out <- list(ACS=ACS,ACSpvalue=ACSpvalue,ADS=ADS,ADSpvalue=ADSpvalue)
  return(out)
}
