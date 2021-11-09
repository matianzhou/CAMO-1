##' The \code{mdsGlobal} is function to generate : multidimensional scaling plot for genome-wide c-scores or d-scores.
##' @title Genome-wide MDS plot for multiple pairs
##' @param acs: pairwise c-scores or d-scores matrix (c-scores is recommended)
##' @param model.name: a vector of dataset names.
##' @param sep: a character string to separate the data terms. Please avoid characters appeared in model.name.
##' @param file: the output file name.
##' @return A MDS plot of c-scores or d-scores for all study pairs to visualize global distance between each data pair.
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step (see the example in function 'merge')
##' #ACS_ADS_global from the multi_ACS_ADS_global step (see the example in 'multi_ACS_ADS_global')
##' dataset.names = c("hb","hs","ht","ha","hi","hl",
##'                   "mb","ms","mt","ma","mi","ml")
##' mdsGlobal(ACS_ADS_global$ACS,dataset.names,sep="_",file="~/globalMDS.pdf")
##' }

mdsGlobal <- function(acs,model.name,sep="_",file) {
  M <- length(model.name)
  acs.value = c(acs)
  names(acs.value) = c(sapply(rownames(acs),function(x) paste(x,rownames(acs),sep=sep)))
  if(M<=2){
    stop("At least three studies are required for Multidimensional Scaling plot...")
  }
  distF <- ACStransform(acs.value,theta=7)
  d <- matrix(NA,nrow=M,ncol=M)
  rownames(d) <- colnames(d) <- model.name

  for(i in 1:(M-1)){
    for(j in (i+1):M){
      name1 <- rownames(d)[i]
      name2 <- rownames(d)[j]
      d[name1,name2] <- d[name2,name1] <- distF[paste(name1,name2,sep=sep)]
    }
  }

  diag(d) <- 0
  dist <- as.dist(d,upper = TRUE, diag = TRUE)
  fit <- sammon(d=dist, y= jitter(cmdscale(dist, 2)), k=2) # k is the number of dim

  x <- fit$points[,1]
  y <- fit$points[,2]
  xlimit <- ifelse(abs(min(x))>abs(max(x)),abs(min(x)),abs(max(x)))
  ylimit <- ifelse(abs(min(y))>abs(max(y)),abs(min(y)),abs(max(y)))

  color <- rainbow(M,s=0.5,v=1,alpha=1)
  pdf(file)
  p<-ggplot() +
    ggtitle("Global MDS plot") +
    xlab("Coordinate 1") + ylab("Coordinate 2") +
    xlim(c(-xlimit-0.5,xlimit+0.5)) + ylim(c(-ylimit-0.5,ylimit+0.5)) +
    geom_point(aes(x, y), color = color  ,size=6) +
    geom_text_repel(aes(x, y, label = rownames(d),fontface="bold"),size=8) +
    theme(text = element_text(face = "bold", color = "black", size=20))
  print(p)
  dev.off()
  return(p)
}
ACStransform <- function(ACS, theta=7) {
  trun.ACS <-ifelse(ACS<0,0,ACS)
  trsf.ACS <- theta*exp(-theta*trun.ACS)
  return(trsf.ACS)
}
