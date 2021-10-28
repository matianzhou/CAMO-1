##' The \code{ACS_ADS_DE_plot} is function to visualize DE evidence and the significance of c-scores and d-scores in one plot. X and y axes represent the average
##' DE posterior probabilities in each study pair, size of dots represent the significance of c-scores (upper right) or d-scores (lower left).
##' @title DE evidence visualization with c-scores and d-scores
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param ACS_ADS_pathway: a list of four data frames: pathway specific c-scores, d-scores
##' and their permuted p-value (rows are pathways and columns are studies).
##' @param dataset.names: a vector of dataset names.
##' @param cluster.lb: a vector of cluster label named by pathway names. optional.
##' @param select.pathway.list: a list of selected pathways (containing gene components).
##' @param highlight.pathway.index: a numeric vector indicating which pathways in the select.pathway.list need to be highlighted.
##' @return A collection of DE evidence plots for each data pair. In each individual panel, each dot is a pathway sized by p(c-score) (orange in upper right)
##' or p(d-score) (blue in lower left).X and y axes are the average absolute DE posterior probabilities.
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step (see the example in function 'merge')
##' #select.pathway.list from the pathSelect step (see the example in function 'pathSelect')
##' #ACS_ADS_pathway from the multi_ACS_ADS_pathway step (see the example in 'multi_ACS_ADS_pathway')
##' #highlight and label the first three pathways in selec.pathway.list:
##' highlight.pathway.index = c(1,2,3)
##' ACS_ADS_DE_plot = function(mcmc.merge.list,ACS_ADS_pathway,dataset.names,
##'                            cluster.lb=NULL,select.pathway.list,
##'                            highlight.pathway.index=highlight.pathway.index,plot.path=NULL)
##' }

ACS_ADS_DE_plot = function(mcmc.merge.list,ACS_ADS_pathway,dataset.names,
                           cluster.lb=NULL,select.pathway.list,
                           highlight.pathway.index=NULL,plot.path=NULL){

  if(is.null(plot.path)){
    plot.path = getwd()
  }

  ACS_pvalue = ACS_ADS_pathway$ACSpvalue.mat
  ACSlog10p.mat = -log10(ACS_pvalue)
  ADS_pvalue = ACS_ADS_pathway$ADSpvalue.mat
  ADSlog10p.mat = -log10(ADS_pvalue)

  P <- nrow(ACSlog10p.mat)
  allgenes <- rownames(mcmc.merge.list[[1]])
  pm.list <- lapply(1:length(mcmc.merge.list), function(x)
    apply(mcmc.merge.list[[x]],1,mean))
  names(pm.list) <- dataset.names


  DEevid <- matrix(NA,P,length(dataset.names))
  rownames(DEevid) = row.names(ACS_pvalue)
  colnames(DEevid) = dataset.names
  for(j in 1:P){
    print(j)
    pathj <- rownames(ACSlog10p.mat)[j]
    genej <- select.pathway.list[[pathj]]
    intergenej <- intersect(genej,allgenes)

    for(ds in dataset.names){
      DEevid[j,ds] <- mean(abs(pm.list[[ds]])[intergenej],na.rm = T)
    }
  }
  write.csv(DEevid,paste0(plot.path,"/DEevid_abs.csv"))


  dpairs = combn(dataset.names,2)
  dpairs.index = combn(length(dataset.names),2)
  Plist = lapply(1:ncol(dpairs), function(x) {
    ds1 <- dpairs[1,x]
    ds2 <- dpairs[2,x]
    DEevid1 = DEevid[,ds1]
    DEevid2 = DEevid[,ds2]
    ACSp <- ACS_pvalue[,paste(ds1,ds2,sep="_")]
    ADSp <- ADS_pvalue[,paste(ds1,ds2,sep="_")]

    #index_pos = index_neg =1:P
    #index_pos[-c(highlight_index_pos,highlight_index_neg)] = ""
    #index_neg[-c(highlight_index_pos,highlight_index_neg)] = ""

    plist = ACS_ADS_DE(ds1,ds2,DEevid1,DEevid2,ACSp,ADSp,cluster.lb,
                       highlight.pathways=highlight.pathway.index,size.scale = 1)
    return(plist)
  })

  Pcordi = cbind(rep(seq(1:length(dataset.names)),each = length(dataset.names)),
                 rep(seq(1:length(dataset.names)),times = length(dataset.names)))
  Plist.org = list()
  for (i in 1:nrow(Pcordi)) {
    print(i)
    if (Pcordi[i,1] < Pcordi[i,2]){
      plist.index = which(sapply(1:ncol(dpairs.index),function(c) {
        all(dpairs.index[,c] == Pcordi[i,])
      }))
      Plist.org[[i]] = Plist[[plist.index]][[1]]
    } else if (Pcordi[i,1] > Pcordi[i,2]){
      plist.index = which(sapply(1:ncol(dpairs.index),function(c) {
        (dpairs.index[1,c] == Pcordi[i,2])&(dpairs.index[2,c] == Pcordi[i,1])
      }))
      Plist.org[[i]] = Plist[[plist.index]][[2]]
    } else {
      Plist.org[[i]] = rectGrob(gp=gpar(fill="white"))
    }
  }
  lay = matrix(seq(1:length(Plist.org)),
               nrow = length(dataset.names),
               ncol = length(dataset.names), byrow = T)

  pdf(paste0(plot.path,"/multiPlot_ACS_ADS_DE_plot.pdf"),width = 50,height = 50)
  grid.arrange(grobs = Plist.org, layout_matrix = lay)
  dev.off()

}
