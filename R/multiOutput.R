##' Downstream visualization tools: visualization outputs including overall
##' pathway clustering and output for each pathway
##' The \code{multiOutput} is function to visualize pathway ACS/ADS results: including overall
##' pathway clustering outputs, comembership heatmaps, model MDS plot, model clustering output,
##' heatmap of gene posterior mean, kegg pathway topology for each pathway.
##' @title Downstream visualization tools
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param dataset.names: a vector of dataset names.
##' @param select.pathway.list: a list of selected pathways (containing gene components).
##' @param ACS_ADS_pathway: a list of four data frames: pathway specific ACS values, ADS values
##' and their permuted p-value (pathway on rows, column being ACS/ADS value or the p-values).
##' @param output: seven options: "clustPathway" (pathway clustering),"comemberPlot" (heatmaps of
##' probabilities that each model pair been clustered together in pathways in each pathway cluster),
##' "mdsModel"(model MDS plot),"clustModel" (model clustering output), "genePM" (heatmap of gene
##' posterior mean), "keggView" (KEGG pathway topology, default is human - hsa, for other species,
##' KEGG organism name and gene Entrez ID needs to be provided as 'KEGG.genes.entrezID'),"reactomeView"
##' (Reactome pathway topology, default is human - HSA, for other species, Reactome organism name needs
##' to be provided as "reactome.species"). In single pair case (2 DTS's), clustering analysis is not
##' applicable and only 'genePM','keggView' can be generated. For details, please refer to manuscript.
##' "output" cannot be empty.
##' @param optK: Optimal number of clusters based on pathway clustering diagnostic results.For
##' "clustPathway" output only.
##' @param sil_cut: silhouette index to control scatterness. Larger value indicates tigher cluster and
##'  more scattered pathways.
##' @param ADS_cluster: whether clustered by ADS in clustPathway. Default is FALSE.
##' @param hashtb: a flat hash table for text mining. Two prepared hastb table are provided: "hashtb_human"
##'  and "hashtb_worm". Please refer to Zeng, Xiangrui, et al. "Comparative Pathway Integrator: a framework
##' of meta-analytic integration of multiple transcriptomic studies for consensual and differential pathway
##' analysis." Genes 11.6 (2020): 696.
##' @param pathways: complete pathway names for text mining.
##' @param keywords_cut: keywords above this cut will be shown in text mining spreadsheet output.
##' @param text.permutation: select from "all" or "enriched". In text mining, "all" permutates pathways from
##' full pathway.list provided while "enriched" permutates from selected pathways. "all" is suitable for
##' cross-species comparision while "enriched" is recommended for within-species comparision.
##' @param ViewPairSelect: which two datasets to view in KEGG/Reactome topology. All pairs will be
##' considered under default (may take a while).
##' @param comemberProb_cut: probability below this cut will be colored blue in comembership heatmaps.
##' @param kegg.species: KEGG species abbreviation. For "keggView" only. Default is "hsa".
##' @param KEGG.genes.entrezID: a data frame which maps gene names in mcmc.merge.list (first column) to
##' Entrez IDs (second column). Only required for non-human genes after merging. For "keggView" only.
##' @param reactome.species: Reactome species abbreviation. For "reactomeView" only. Default is "HSA".
##' @param Reactome.genes.topologyGtype: a data frame which maps gene names in mcmc.merge.list (first
##' column) to gene name types in Reactome topology (second column). For "reactomeView" only.
##'
##' @return stored output in created folders.
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step (see the example in function 'merge')
##' #select.pathway.list from the pathSelect step (see the example in function 'pathSelect')
##' #ACS_ADS_pathway from the multi_ACS_ADS_pathway step (see example in 'multi_ACS_ADS_pathway')
##' data(hashtb_human) #include hashtb & pathways for text mining
##' dataset.names = c("hb","hs","ht","ha","hi","hl",
##'                   "mb","ms","mt","ma","mi","ml")
##' #1. step1: select K by elbow plot from consensus clustering
##' ACSpvalue.mat = ACS_ADS_pathway$ACSpvalue.mat
##' results = clustPathway(ACSpvalue.mat)
##'
##' #2. step2: run multiOutput with pre-selected K
##' multiOutput(mcmc.merge.list,dataset.names,select.pathway.list,ACS_ADS_pathway,
##'            output=c("clustPathway","comemberList","mdsModel","clustModel","genePM","keggView"),
##'            hashtb=hashtb,pathways=pathways,optK = K,thres = 0.2)
##' }
multiOutput <- function(mcmc.merge.list,dataset.names,select.pathway.list,ACS_ADS_pathway,
                        output=c("clustPathway","mdsModel","clustModel","genePM","keggView","reactomeView",
                                 "comemberPlot"),
                        optK=NULL,sil_cut=0.1,ADS_cluster=FALSE,hashtb=NULL,pathways=NULL,keywords_cut=0.05,
                        text.permutation = "all",comemberProb_cut=0.7,
                        ViewPairSelect = NULL,kegg.species="hsa",KEGG.genes.entrezID=NULL,
                        reactome.species="HSA",Reactome.genes.topologyGtype=NULL) {#ViewPairSelect:a subset of dataset.names for keggView

  ### Multiple pairs output including the following outputs:
  ## pathway clustering (including consensus clust, heatmap, mds, pathway text mining)
  ## & pathway level: per pathway mdsModel plot, clustModel heatmap, gene heatmap,
  ## pathview (KEGG only)

  if(length(output)==0 || is.null(output)){
    stop("at least one type of output has to be chosen")
  }

  if(("comemberPlot" %in% output) & (!"clustPathway" %in% output)){
    print("comemberPlot is depend on pathway clustering results, can not generate without 'clustPathway'...")
  }

  orig.path <- getwd()
  pathway.name <- names(select.pathway.list)
  K <- length(pathway.name)

  if(ADS_cluster == TRUE){
    AS.mat = ACS_ADS_pathway$ADS.mat
    ASpvalue.mat = ACS_ADS_pathway$ADSpvalue.mat
  }else{
    AS.mat = ACS_ADS_pathway$ACS.mat
    ASpvalue.mat = ACS_ADS_pathway$ACSpvalue.mat
  }


  M <- length(mcmc.merge.list)
  P <- ncol(AS.mat)

  if( (M<=2) & (length(intersect(output,c("clustPathway","mdsModel","clustModel","comemberPlot"))) != 0)){
    print("MDS and clustering are not applicable to 2 DTS's case, removed from 'output'...")
    output = intersect(output,c('genePM','keggView'))
  }

  if("clustPathway" %in% output){
    dir.path <- "clustPathway"
    if (!file.exists(dir.path)) dir.create(dir.path)
    setwd(paste(orig.path,"/",dir.path,sep=""))

    #1. consensus clustering
    results <- clustPathway(ASpvalue.mat)

    #2. determine optimal number of clusters
    if(is.null(optK)){
      optK <- clustNumber(results)
    }
    cluster.assign <- results[[optK]]$consensusClass

    #3. identify scattered objects
    scatter.index <- scatter(-log10(ASpvalue.mat), cluster.assign, sil_cut=sil_cut)
    cluster = list(cluster.assign=cluster.assign,scatter.index=scatter.index)
    save(cluster,file = "cluster_labels.RData")

    #4. mds plot
    res <- mdsPathway(acsPvalue=ASpvalue.mat,
                      cluster.assign=cluster.assign,
                      scatter.index=scatter.index)
    #5. text mining
    if(!is.null(hashtb) && !is.null(pathways)){
      textMine(hashtb=hashtb,pathways=pathways,cluster.assign=cluster.assign,scatter.index=scatter.index,thres=keywords_cut,permutation=text.permutation)
    }

    #6. heatmap
    res <- heatmapPathway(acsPvalue=ASpvalue.mat,
                          cluster.assign=cluster.assign,
                          scatter.index=scatter.index)

    #7. ACS_ADS_DE plot (colored by cluster)
    res <- ACS_ADS_DE_cluster(mcmc.merge.list=mcmc.merge.list,ACS_ADS_pathway=ACS_ADS_pathway,
                              dataset.names=dataset.names,select.pathway.list=select.pathway.list,
                              cluster.assign=cluster.assign,scatter.index=scatter.index,plot.path=NULL)


    if("comemberPlot" %in% output){
      print("Construct comembership matrix...")
      dir.path <- "comemberPlot"
      if (!file.exists(paste(orig.path,"/",dir.path,sep=""))) dir.create(paste(orig.path,"/",dir.path,sep=""))
      setwd(paste(orig.path,"/",dir.path,sep=""))

      model.cluster.result <- vector("list",length=K)
      names(model.cluster.result) <- rownames(ASpvalue.mat)

      for(k in 1:K) {
        model.cluster.result[[k]] <- SA_algo(unlist(c(ASpvalue.mat[k,])),dataset.names,sep="_")
      }
      pathway.cluster.assign = cluster.assign
      pathway.cluster.assign[scatter.index] = "scatter"

      if(is.null(scatter.index)){
        Cvec = 1:optK
      }else{
        Cvec = c(1:optK,"scatter")
      }

      comember.list <- vector("list",length=length(Cvec))
      for (c in 1:length(Cvec)){
        select.pathways <- names(pathway.cluster.assign)[pathway.cluster.assign==Cvec[c]]
        denom <- length(select.pathways)
        model.result.pathways <- matrix(unlist(model.cluster.result[select.pathways]),nrow=denom,ncol=M,byrow =T )
        rownames(model.result.pathways) <- select.pathways
        colnames(model.result.pathways) <- dataset.names
        comember.mat <- matrix(1,M,M)
        rownames(comember.mat) <- colnames(comember.mat) <- dataset.names
        for (i in 1:(M-1)){
          for (j in (i+1):M){
            if (denom==1){
              comember.mat[j,i] = 1
            }else{
              model1 <- dataset.names[i]
              model2 <- dataset.names[j]
              twomodel.result <- model.result.pathways[,c(model1,model2)]
              comember.mat[j,i] <- comember.mat[i,j] <- sum(apply(twomodel.result,1,function(x) x[1]==x[2]))/denom
            }
          }
        }
        comember.list[[c]] <- comember.mat
      }
      names(comember.list) = Cvec
      save(comember.list,file="comember.list.RData")

      #1. thresholding:
      threshold = comemberProb_cut

      mat.thres <- function(mat, threshold){
        mat[which(mat <= threshold)] <- 0
        return(mat)
      }

      #2. matrix multiplication


      source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

      for(i in 1:length(comember.list)){
        mat <- comember.list[[i]]

        par(font.main=2)
        par(font.axis=2)
        par(xpd=FALSE)

        #thresholding

        thres.mat <- mat.thres(mat,threshold)

        pdf(paste("ComemMat_cluster_",names(comember.list)[i],"_threshold_",threshold,".pdf",sep=""))


        hm <- heatmap.3(thres.mat,
                        main=NULL,
                        cexCol=1.2,cexRow=1.2,
                        colsep=1:nrow(mat),
                        rowsep=1:ncol(mat),
                        sepwidth=c(0.02, 0.02),  # width of the borders
                        sepcolor='black',
                        symbreaks=T,key=T, keysize=1,symkey=F,
                        dendrogram=c('none'),density.info="none",
                        trace="none",Rowv=T,Colv=T,symm=F,
                        #srtCol=50,
                        col=bluered,breaks=seq(0,1,by=0.01),
                        margins=c(12,14))

        dev.off()
      }
    }
  }
  setwd(orig.path)

  if("mdsModel" %in% output){
    dir.path <- "mdsModel"
    if (!file.exists(dir.path)) dir.create(dir.path)
    setwd(paste(orig.path,"/",dir.path,sep=""))

    for(k in 1:K){
      print(paste("mdsModel",k,sep=":"))
      pathk.name <- pathway.name[k]
      res <- mdsModel(unlist(c(AS.mat[k,])),dataset.names,pathk.name,sep="_")
    }
  }
  setwd(orig.path)

  if("clustModel" %in% output){
    dir.path <- "clustModel"
    if (!file.exists(dir.path)) dir.create(dir.path)
    setwd(paste(orig.path,"/",dir.path,sep=""))

    for(k in 1:K){
      print(paste("clustModel",k,sep=":"))
      pathk.name <- pathway.name[k]
      cluster.assign.path <- try(SA_algo(unlist(c(ASpvalue.mat[k,])),dataset.names,sep="_"))
      if(length(unique(cluster.assign.path))>1 && class(cluster.assign.path) != "try-error" ){
        res <- clustModel(unlist(c(ASpvalue.mat[k,])),dataset.names,cluster.assign.path,
                          pathk.name,sep="_")
      }
      if(length(unique(cluster.assign.path))==1 || class(cluster.assign.path) == "try-error" ){
        warning("clustModel only identifies one cluster")
        res <- clustModelOne(unlist(c(ASpvalue.mat[k,])),dataset.names,
                             pathk.name,sep="_")
      }
    }
  }
  setwd(orig.path)


  if("genePM" %in% output){
    dir.path <- "genePM"
    if (!file.exists(dir.path)) dir.create(dir.path)
    setwd(paste(orig.path,"/",dir.path,sep=""))
    for(k in 1:K){
      print(paste("genePM",k,sep=":"))
      pathk.name <- pathway.name[k]
      pathway.genes <- select.pathway.list[[k]]
      signPM.list <- lapply(mcmc.merge.list,function(x) apply(x,1,mean))
      names(signPM.list) <- dataset.names
      hm <- genePM(signPM.list, pathway.genes=pathway.genes,
                   pathway.name=pathk.name)
    }
  }
  setwd(orig.path)

  if("keggView" %in% output){
    if(sum(grepl("KEGG",pathway.name))==0) {
      warning("No KEGG pathways")
    } else{
      dir.path <- "keggView"
      if (!file.exists(dir.path)) dir.create(dir.path)
      setwd(paste(orig.path,"/",dir.path,sep=""))

      map.ls = as.list(as.list(org.Hs.egALIAS2EG))

      kegg_id_name = lapply(keggList("pathway",kegg.species),function(x) strsplit(x," - ")[[1]][-length(strsplit(x," - ")[[1]])])

      if(is.null(ViewPairSelect)){
        data.pair = combn(dataset.names,2)
      }else{
        data.pair = combn(ViewPairSelect,2)
      }
      select.kegg.path.name = pathway.name[grep("KEGG",pathway.name)]## KEGG pathways selected
      K_KEGG = length(select.kegg.path.name)
      for(k in 1:K_KEGG){
        print(paste("keggView",k,sep=":"))
        keggk.name = select.kegg.path.name[k]
        for (i in 1:ncol(data.pair)) {
          dat1.name = data.pair[1,i]
          dat2.name = data.pair[2,i]
          dat1 = mcmc.merge.list[[dat1.name]]
          dat2 = mcmc.merge.list[[dat2.name]]
          overlap.genes <- intersect(rownames(dat1),select.pathway.list[[keggk.name]])
          signPM.mat <- cbind(apply(dat1[overlap.genes,],1,mean),
                              apply(dat2[overlap.genes,],1,mean))
          colnames(signPM.mat) <- c(dat1.name,dat2.name)


          std.genes <- rownames(signPM.mat)
          if(is.null(KEGG.genes.entrezID)){
            rownames(signPM.mat) = sapply(std.genes,function(g) map.ls[[g]][[1]])
          }else{
            entrezID = KEGG.genes.entrezID[match(std.genes,KEGG.genes.entrezID[,1]),2]
            na.index = which(is.na(entrezID))
            if(length(na.index != 0)){
              entrezID = entrezID[-na.index]
              signPM.mat = signPM.mat[-na.index,]
            }
            rownames(signPM.mat) = entrezID
          }

          keggk.name.clean = gsub("KEGG ","",keggk.name)
          pathwayID = gsub(paste0("path:",kegg.species),"",names(kegg_id_name)[which(kegg_id_name == keggk.name.clean)])

          res = pathview(gene.data = signPM.mat, pathway.id = pathwayID,
                         species = kegg.species, out.suffix = "", kegg.native = T,
                         key.pos = "bottomright", map.null=T,cex = 0.15)


          keggID = paste(kegg.species,pathwayID,sep="")
          keggk.name0 = gsub(" / ","_",keggk.name,fixed = T)
          file.rename(paste(keggID,"..multi.png",sep=""),
                      paste(keggk.name0,"_",dat1.name,"_",dat2.name,".png",sep=""))
          file.remove(paste(keggID,".xml",sep=""))
          file.remove(paste(keggID,".png",sep=""))
        }
      }
    }
  }
  setwd(orig.path)

  if("reactomeView" %in% output){
    if(sum(grepl("Reactome",pathway.name))==0) {
      warning("No Reactome pathways")
    } else{
      dir.path <- "reactomeView"
      if (!file.exists(dir.path)) dir.create(dir.path)
      setwd(paste(orig.path,"/",dir.path,sep=""))

      source_python(system.file("ImageProcess.py", package = "CAMO"))
      file.copy(from = system.file("pallete.jpeg", package = "CAMO"),
                to   = getwd())

      #genes in mcmc.merge.list, pathway.list should both be the type shown on Reactome topology
      if(!is.null(Reactome.genes.topologyGtype)){
        #match mcmc.merge.list genes
        mcmcG = row.names(mcmc.merge.list[[1]])
        mcmcG.topology = Reactome.genes.topologyGtype[match(mcmcG,Reactome.genes.topologyGtype[,1]),2]
        na.index = which(is.na(mcmcG.topology))
        if(length(na.index) !=0){
          mcmcG.topology.rmna = mcmcG.topology[-na.index]
          mcmc.merge.list.topology = lapply(mcmc.merge.list, function(x) {
            x_rm = x[-na.index,]
            row.names(x_rm) = mcmcG.topology.rmna
            return(x_rm)
          })
        }else{
          mcmc.merge.list.topology = lapply(mcmc.merge.list, function(x) {
            row.names(x) = mcmcG.topology
            return(x)
          })
        }
        #match pathway genes
        select.pathway.list.topology = lapply(select.pathway.list, function(x){
          pathwayG.topology = Reactome.genes.topologyGtype[match(x,Reactome.genes.topologyGtype[,1]),2]
          pathwayG.topology = pathwayG.topology[-which(is.na(pathwayG.topology)|pathwayG.topology == "")]
          return(pathwayG.topology)
        })
      }else{
        mcmc.merge.list.topology = mcmc.merge.list
        select.pathway.list.topology = select.pathway.list
      }

      pathid2name = as.list(reactomePATHID2NAME)
      pathid2name_species = pathid2name[grep(reactome.species,names(pathid2name))]
      pathid2name_species_clean = lapply(pathid2name_species, function(x) strsplit(x,": ")[[1]][2])

      pathway.name.topology = names(select.pathway.list.topology)
      select.path.name = gsub("Reactome ","",pathway.name.topology[grep("Reactome",pathway.name.topology)])
      select.path.name.intersect = intersect(pathid2name_species_clean,select.path.name)
      select.path.id.intersect = names(pathid2name_species_clean)[match(select.path.name.intersect,pathid2name_species_clean)]
      print(paste0("Matched ",length(select.path.id.intersect)," Reactome pathway IDs"))

      select.path.name.intersect2 = paste0("Reactome ",select.path.name.intersect)
      reactome.index = match(select.path.name.intersect2, pathway.name.topology)
      reactome.df = data.frame(ID = as.character(select.path.id.intersect),
                               Genes = as.character(sapply(select.pathway.list.topology[select.path.name.intersect2], function(x) paste(x,collapse = " "))))

      if(is.null(ViewPairSelect)){
        data.pair = combn(dataset.names,2)
      }else{
        data.pair = combn(ViewPairSelect,2)
      }
      for(i in 1:length(select.path.id.intersect)){
        print(paste("reactomeView",i,sep=":"))
        aID = select.path.id.intersect[i]
        aname = select.path.name.intersect2[i]
        if(grepl(":|/", aname)) {
          aname1 = gsub(":|/", "-", aname)
        } else {
          aname1 = aname
        }
        dir.path <- paste0("reactomeView/", aname1)
        if (!file.exists(dir.path)) dir.create(dir.path)
        setwd(paste(orig.path,"/",dir.path,sep=""))
        file.copy(from = system.file("pallete.jpeg", package = "CAMO"),
                  to   = getwd())
        for (j in 1:ncol(data.pair)) {
          tryCatch({
            #aID = select.path.id.intersect[i]
            #aname = select.path.name.intersect2[i]

            dat1.name = data.pair[1,j]
            dat2.name = data.pair[2,j]
            dat1 = mcmc.merge.list.topology[[dat1.name]]
            dat2 = mcmc.merge.list.topology[[dat2.name]]
            print(paste0(aID," ",dat1.name," ",dat2.name))
            overlap.genes0 = intersect(rownames(dat1),select.pathway.list.topology[[aname]])
            if(length(overlap.genes0) != 0){
              signPM.mat = cbind(apply(dat1[overlap.genes0,],1,mean),
                                 apply(dat2[overlap.genes0,],1,mean))
              signPM = data.frame(Genes = row.names(signPM.mat),signPM.mat)
              colnames(signPM) = c("Genes","dat1","dat2")


              CallFromR(signPM=signPM,ReactomePath=reactome.df,datadir=paste0(getwd(),"/"),pathwayID=aID)

              file.rename(paste(aID,"New.jpeg",sep=""),
                          paste(aID,"_",dat1.name,"_",dat2.name,".jpeg",sep=""))
              file.remove(paste(aID,".jpeg",sep=""))
              file.remove(paste(aID,"(1).jpeg",sep=""))
              file.remove(paste(aID,".sbgn",sep=""))
            }
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
        file.remove("pallete.jpeg")
        setwd(orig.path)

      }
    }
  }
  setwd(orig.path)

  print("Multiple pairs analysis completed.")

}

#internal functions for multiOutput:

clustPathway <- function(acsPvalue) {
  #require(ConsensusClusterPlus)
  #set your working dir, automatically save there

  results = ConsensusClusterPlus(d=t(-log10(acsPvalue)),maxK=10,reps=50,pItem=0.8,pFeature=1,title="Pathway clustering",clusterAlg="hc",innerLinkage="ward.D2",finalLinkage="ward.D2",seed=15213,plot="png")

  return(results) ## a list of K elements (each k represents number of clusters)
}

clustNumber <- function(results){
  Kvec = 2:length(results);
  x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC = rep(NA,length(Kvec))
  names(PAC) = paste("K=",Kvec,sep="") # from 2 to 10 (maxK)
  for(i in Kvec){
    M = results[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }#end for i
  # The optimal K
  optK = Kvec[which.min(PAC)-2]
  return(optK)
}

ACS_ADS_DE_cluster = function(mcmc.merge.list,ACS_ADS_pathway,dataset.names,
                              select.pathway.list,cluster.assign,scatter.index,
                              plot.path=NULL){

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
    pathj <- rownames(ACSlog10p.mat)[j]
    genej <- select.pathway.list[[pathj]]
    intergenej <- intersect(genej,allgenes)

    for(ds in dataset.names){
      DEevid[j,ds] <-  mean(abs(pm.list[[ds]])[intergenej],na.rm=T)
    }
  }
  write.csv(DEevid,paste0(plot.path,"/DEevid_abs.csv"))

  cluster.lb = cluster.assign
  cluster.lb[scatter.index] = "scatter"
  K = length(unique(cluster.assign))

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

    plist = ACS_ADS_DE(ds1,ds2,DEevid1,DEevid2,ACSp,ADSp,cluster.lb,size.scale = 1)
    return(plist)
  })

  Pcordi = cbind(rep(seq(1:length(dataset.names)),each = length(dataset.names)),
                 rep(seq(1:length(dataset.names)),times = length(dataset.names)))
  Plist.org = list()
  for (i in 1:nrow(Pcordi)) {
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

  pdf(paste0(plot.path,"/multiPlot_ACS_ADS_DE_K",K,".pdf"),width = 50,height = 50)
  grid.arrange(grobs = Plist.org, layout_matrix = lay)
  dev.off()

}

mdsPathway <- function(acsPvalue,cluster.assign,scatter.index=NULL) {

  ## plot MDS for all pathways
  ## acsPvalue is a matrix of K (pathways) rows and choose(M,2) columns
  ## cluster.assign is the result from consensus clustering

  C <- length(unique(cluster.assign)) #number of clusters
  if(C>9){
    warning("Too many clusters, not enough colors")
  }

  dist.mat <-  dist(-log10(acsPvalue),method = "euclidean",
                    upper = TRUE, diag = TRUE)

  fit <- cmdscale(dist.mat,k=2)
  x <- fit[,1]
  y <- fit[,2]
  xlimit <- ifelse(abs(min(x))>abs(max(x)),abs(min(x)),abs(max(x)))
  ylimit <- ifelse(abs(min(y))>abs(max(y)),abs(min(y)),abs(max(y)))
  xcenter <- tapply(x,as.factor(cluster.assign),mean)
  ycenter <- tapply(y,as.factor(cluster.assign),mean)

  if(!is.null(scatter.index)){
    cluster.assign[scatter.index] <- -1
  }

  unique.color <- rainbow(C)
  unique.shape <- 1:C
  sizes <- shapes <- colors <- cluster.assign
  for(i in 1:(C+1)){
    colors[cluster.assign==i] <- unique.color[i]
    shapes[cluster.assign==i] <- unique.shape[i]
    sizes[cluster.assign==i] <- 2
    if(i== (C+1)){
      colors[cluster.assign== -1] <- "gray50"
      shapes[cluster.assign== -1] <- 20
      sizes[cluster.assign== -1] <- 1
    }
  }
  #pdf(paste("mdsPathway","_K_",C,".pdf",sep=""))
  jpeg(paste("mdsPathway","_K_",C,".jpeg",sep=""),quality = 100)

  p <- ggplot() +
    ggtitle("") +
    xlab("Coordinate 1") + ylab("Coordinate 2") +
    xlim(c(-xlimit,xlimit)) + ylim(c(-ylimit,ylimit)) +
    geom_point(aes(x, y), shape=shapes,
               color = colors ,size=sizes) +
    geom_point(aes(xcenter,ycenter),
               shape=unique.shape, color = unique.color,
               size =5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 15, hjust=0.5,face="bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  print(p)
  dev.off()
  return(p)
}

heatmapPathway <- function(acsPvalue, cluster.assign,scatter.index=NULL){

  ## cluster.assign is the result from consensus clustering
  acsPvalue <- data.matrix(acsPvalue)
  C <- length(unique(cluster.assign)) #number of clusters
  dataOrder <- -log10(acsPvalue)[unlist(sapply(1:C,function(x) which(cluster.assign==x))),]
  colnames(dataOrder) <- colnames(acsPvalue)
  row.sep <-  c(0,cumsum(unlist(sapply(1:C,function(x) sum(cluster.assign==x)))) )

  if(!is.null(scatter.index)){
    cluster.assign[scatter.index] <- -1
    dataOrder <- -log10(acsPvalue)[unlist(sapply(c(1:C,-1),function(x) which(cluster.assign==x))),]
    colnames(dataOrder) <- colnames(acsPvalue)
    row.sep <-  c(0,cumsum(unlist(sapply(c(1:C,-1),function(x) sum(cluster.assign==x)))) )
  }
  #rownames(dataOrder) <- sapply(rownames(dataOrder),function(x) substr(x,1,15))
  ordered.cluster.assign <- cluster.assign[rownames(dataOrder)]
  row.colors <- rep(NA, length(ordered.cluster.assign) )
  for(i in 1:length(row.colors)){
    if(ordered.cluster.assign[i]== -1){
      row.colors[i] <- "gray"
    } else {
      row.colors[i] <- rainbow(C)[ordered.cluster.assign[i]]
    }
  }
  #pdf(paste("heatmapPathway","_K_",C,".pdf",sep=""))
  jpeg(paste("heatmapPathway","_K_",C,".jpeg",sep=""),quality = 100)
  par(cex.main=1, font.lab=2, font.axis=2)
  hm<-heatmap.2(dataOrder, symm=F,main=NULL,
                cexCol=0.7,cexRow=0.3,adjCol= c(NA,-1),
                rowsep=row.sep,
                sepwidth=c(0.1, 0.3),  # width of the borders
                sepcolor=c('white'),scale='none',
                symbreaks=T,key=T, keysize=1,symkey=F,
                dendrogram=c('none'),density.info="none",
                trace="none",Rowv=F,Colv=T,
                srtCol=50,RowSideColors=row.colors,
                col=greenred,breaks=seq(0,max(dataOrder),by=0.01),
                key.ytickfun=function(){
                  side = 2
                } )
  dev.off()

  return(hm)

}

SA_algo <- function(acsPvaluePath,model.name,sep,Tm=10,P=0.5,
                    mu=0.9,epsilon=1e-5,N=1000,seed=12345){
  ## Total possible configurations = choose(n,1) + choose(n,2) + ...
  ## clustering models using SA algorithm
  ## delta.mat is a matrix of pairwise -log10(acsPvalue) of M rows and M columns
  ## with model names
  M <- length(model.name)
  distF <- -log10(acsPvaluePath)
  delta.mat <- matrix(NA,nrow=M,ncol=M)
  rownames(delta.mat) <- colnames(delta.mat) <- model.name

  for(i in 1:(M-1)){
    for(j in (i+1):M){
      name1 <- rownames(delta.mat)[i]
      name2 <- rownames(delta.mat)[j]
      delta.mat[name1,name2] <- delta.mat[name2,name1] <- distF[paste(name1,name2,sep=sep)]
    }
  }
  diag(delta.mat) <- max(delta.mat,na.rm=T)

  n <- nrow(delta.mat) # =M
  obs.name <- rownames(delta.mat)

  ## initialize
  hc <- hclust(d=dist((max(delta.mat)-delta.mat)^2))
  K <- 3
  a <- cutree(hc, k = K)
  names(a) <- obs.name
  delta.est <- Est_mean(delta.mat,a)
  names(delta.est) <- c(0,1:K)

  count <- 0 #initial
  Jc <- E_tot(delta.mat,a,delta.est)
  pi <- exp(-Jc/Tm) ## Boltzmann dist

  while((count < N) && (Tm >= epsilon)) {
    ##New trial
    u <- runif(1)
    if(u>0.5){
      a_new <- Split(a)
    } else{
      a_new <- Relocate(a)
    }
    names(a_new) <- obs.name
    K_new <- length(unique(a_new))
    delta.est_new <- Est_mean(delta.mat,a_new)
    Jn <- E_tot(delta.mat, a_new, delta.est_new)
    pi_new <- exp(-Jn/Tm)

    if(Jn < Jc) {
      ##accept
      Jc <- Jn;
      a <- a_new;
      K <- K_new;
      delta.est <- delta.est_new;
    } else {
      count <- count + 1;
      r <- min(1,pi_new/pi); ## acceptance prob.
      u <- runif(1);
      if(u>r) {
        ##not accept
        Tm <- Tm*mu
      } else {
        Jc <- Jn;
        a <- a_new;
        K <- K_new;
        delta.est <- delta.est_new;
      }
    }
  }
  return(a) #cluster assigned
}

clustModel <- function(acsPvaluePath,model.name, cluster.assign,pathway.name,sep){
  ## delta.mat is a matrix of pairwise -log(acsPvalue) of M rows and M columns
  ## with model names
  ## cluster.assign: results from SA algorithm
  ## pathway.name: pathway name

  M <- length(model.name)
  distF <- -log10(acsPvaluePath)
  delta.mat <- matrix(NA,nrow=M,ncol=M)
  rownames(delta.mat) <- colnames(delta.mat) <- model.name

  for(i in 1:(M-1)){
    for(j in (i+1):M){
      name1 <- rownames(delta.mat)[i]
      name2 <- rownames(delta.mat)[j]
      delta.mat[name1,name2] <- delta.mat[name2,name1] <- distF[paste(name1,name2,sep=sep)]
    }
  }
  diag(delta.mat) <- max(delta.mat,na.rm=T)

  if(grepl("/",pathway.name)){
    pathway.name <- gsub("/","-",pathway.name)
  }
  if(max(delta.mat)<1){
    breaks=seq(0,max(delta.mat),length.out = 10)
  }else{
    breaks=seq(0,round(max(delta.mat)),by=0.01)
  }
  png(paste(pathway.name,'.png',sep="_"))
  hm <- heatmap.2(delta.mat[order(cluster.assign),order(cluster.assign)],
                  main=pathway.name,
                  cexCol=1,cexRow=1,
                  colsep=cumsum(table(cluster.assign)),
                  rowsep=cumsum(table(cluster.assign)),
                  sepwidth=c(0.05, 0.05),  # width of the borders
                  sepcolor=c('white'),
                  symbreaks=T,key=T, keysize=1,symkey=F,
                  dendrogram=c('none'),density.info="none",
                  trace="none",Rowv=F,Colv=F,
                  srtCol=50, symm=F,
                  col=greenred,breaks = breaks)
  dev.off()
  return(hm)
}

clustModelOne <- function(acsPvaluePath,model.name, pathway.name,sep){
  ## delta.mat is a matrix of pairwise -log(acsPvalue) of M rows and M columns
  ## with model names
  ## cluster.assign: results from SA algorithm
  ## pathway.name: pathway name

  M <- length(model.name)
  distF <- -log10(acsPvaluePath)
  delta.mat <- matrix(NA,nrow=M,ncol=M)
  rownames(delta.mat) <- colnames(delta.mat) <- model.name

  for(i in 1:(M-1)){
    for(j in (i+1):M){
      name1 <- rownames(delta.mat)[i]
      name2 <- rownames(delta.mat)[j]
      delta.mat[name1,name2] <- delta.mat[name2,name1] <- distF[paste(name1,name2,sep=sep)]
    }
  }
  diag(delta.mat) <- max(delta.mat,na.rm=T)
  if(max(delta.mat)<1){
    breaks=seq(0,max(delta.mat),length.out = 10)
  }else{
    breaks=seq(0,round(max(delta.mat)),by=0.01)
  }
  if(grepl("/",pathway.name)){
    pathway.name <- gsub("/","-",pathway.name)
  }

  png(paste(pathway.name,'.png',sep="_"))
  hm <- heatmap.2(delta.mat,
                  main=pathway.name,
                  cexCol=1,cexRow=1,
                  symbreaks=T,key=T, keysize=1,symkey=F,
                  dendrogram=c('none'),density.info="none",
                  trace="none",Rowv=F,Colv=F,
                  srtCol=50, symm=F,
                  col=greenred,breaks=seq(0,round(max(delta.mat)),by=0.01) )
  dev.off()
  return(hm)
}

mdsModel <- function(acsPath,model.name,pathway.name,sep) {
  ## for each pathway, plot MDS for all models
  ## acsP is a vector of choose(M,2) elements (named in paste(name1,name2,sep=""))
  ## model.name is a vector of model names
  M <- length(model.name)
  distF <- ACStransform(acsPath)
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

  if(grepl("/",pathway.name)){
    pathway.name <- gsub("/","-",pathway.name)
  }

  color <- rainbow(M,s=0.5,v=1,alpha=1)
  png(paste(pathway.name,".png",sep=""))
  p<-ggplot() +
    ggtitle(pathway.name) +
    xlab("Coordinate 1") + ylab("Coordinate 2") +
    xlim(c(-xlimit-0.5,xlimit+0.5)) + ylim(c(-ylimit-0.5,ylimit+0.5)) +
    geom_point(aes(x, y), color = color  ,size=6) +
    geom_text_repel(aes(x, y, label = rownames(d),fontface="bold"),size=8) +
    theme(plot.title = element_text(size = 15, hjust=0.5,face="bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  print(p)
  dev.off()
}

genePM <- function(signPM.list, pathway.genes, pathway.name){

  M <- length(signPM.list)
  model.name <- names(signPM.list)
  std.genes <- intersect(names(signPM.list[[1]]),pathway.genes)
  G <- length(std.genes)

  mat <- matrix(0,nrow=G,ncol=M)
  rownames(mat) <- std.genes
  colnames(mat) <- model.name

  for (m in 1:M){
    mat[std.genes,m] <- signPM.list[[m]][std.genes]
  }

  if(grepl("/",pathway.name)){
    pathway.name <- gsub("/","-",pathway.name)
  }

  #pdf(paste(pathway.name,'.pdf',sep=""))
  jpeg(paste(pathway.name,".jpeg",sep=""),quality = 100)
  par(cex.main=1)
  hm <- heatmap.2(mat, symm=F,main=pathway.name,
                  cexCol=1,cexRow=0.6,
                  colsep=c(1:ncol(mat)),
                  rowsep=c(1:nrow(mat)),
                  sepwidth=c(0.01, 0.2),  # width of the borders
                  sepcolor=c('black'),scale='none',
                  symbreaks=T,key=T, keysize=1,symkey=F,
                  dendrogram=c('row'),density.info="none",
                  trace="none",Rowv=T,Colv=F,
                  col=bluered,breaks=seq(-1,1,by=0.001),
                  srtCol=0,adjCol = c(NA,0.5))
  dev.off()
  return(hm)
}

ACStransform <- function(ACS, theta=7) {
  trun.ACS <-ifelse(ACS<0,0,ACS)
  trsf.ACS <- theta*exp(-theta*trun.ACS)
  return(trsf.ACS)
}
