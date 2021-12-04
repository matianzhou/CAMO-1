##' The \code{KEGG_module} is a function to identify concordant or discordant subnetworks in KEGG pathways
##' based on topological regulatory information by generating local module memberships and
##' corresponding p-values at different module sizes.
##' @title KEGG local concordant/discordant module detection algorithm: local module memberships and
##' corresponding p-values at different module sizes
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param dataset.names: a vector of dataset names matched with the mcmc.merge.list.
##' @param KEGGspecies: the KEGG species abbreviation. Default is "hsa".
##' @param KEGGpathwayID: a KEGG pathway ID, not including the organism prefix.
##' @param KEGG.dataGisTopologyG: whether gene names in data are same as entries on KEGG topology. If TRUE,
##' search topology nodes/entries by data gene names directly. Default is FALSE.
##' @param KEGG.dataG2topologyG: a data frame which maps gene names in mcmc.merge.list (first column) to
##' entries on KEGG topology (second column). If NULL & KEGG.dataGisEntrezID=F & KEGGspecies is "hsa", "mmu" or
##' "rno", gene symbols will be automatically mapped to EntrezID by Bioconductor packages
##' "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db". If NULL & KEGG.dataGisEntrezID=F & KEGGspecies="cel",
##' gene symbols will be automatically mapped to WormBase sequence name by Bioconductor package "biomaRt" with prefix "CELE_"
##' added to match entry names in KEGG topology for Caenorhabditis elegans. If NULL & KEGG.dataGisEntrezID=F & KEGGspecies="dme",
##' gene symbols will be automatically mapped to EntrezID and then FlyBase CG IDs by Bioconductor package "org.Dm.eg.db" with
##' prefix "Dmel_" added to match entry names in KEGG topology for Drosophila melanogaster.
##' @param data.pair: a character vector of two study names.
##' @param gene_type: the type of module of interests. This should be one of "concordant" or
##' "discordant".
##' @param DE_PM_cut: only concordant/discordant genes with posterior mean of DE indicators above
##' this value will be considered when searching for modules.
##' @param minM: the miminum module size to consider during searching.
##' @param maxM: the maximum module size to consider during searching. If NULL, maximum module size will be the number of all
##' concordant/discordant genes.
##' @param B: the number of permutations.
##' @param cores: the number of cores to use in permutation (mc.cores parameter in 'mclapply'
##' function).
##' @param search_method: the method used to search modules with small average shortest path. This
##' should be one of "Exhaustive" or "SA" (Simulated-Annealing).
##' @param reps_eachM: the number of searching repetitions at each module size when SA is selected.
##' @param topG_from_previous: the number of top module results stored as initials for next module
##' size when SA is selected.
##' @param Tm0: SA parameter - the initial temparature.
##' @param mu: SA parameter - the temparature multiplier.
##' @param epsilon: SA parameter - the final temparature.
##' @param N: SA parameter - the number of maximum annealing times.
##' @param Elbow_plot: a logical value indicating if an elbow plot of -log10(p-value) at each module
##' size will be saved.
##' @param filePath: the path to save the elbow plot. Default is the current working directory.
##' @param seed: permutation seed.
##' @return A list containing 5 elements.
##' \itemize{
##' \item minG.ls: contains the following information for each module size from minM to maxM.
##' \code{minG} has genes in the module whose average shortest path is optimized. \code{p.mean},
##' \code{p.sd} and \code{sp} are p-values, corresponding standard deviation (sd) and the average
##' shortest path respectively. \code{null.sp.mean} and \code{null.sp.median} are from permutated
##' null distribution. If \code{SA}is selected, the top \code{topG_from_previous} results at each
##' module size is stored in \code{top.G}.
##' \item bestSize: minG.ls results for the largest module size within 2 sd of the smallest p-value.
##' \item mergePMmat: a merged posterior DE mean matrix based on topology nodes (one node can have
##' multiple genes).
##' \item KEGGspecies: the KEGG species abbreviation.
##' \item KEGGpathwayID: the KEGG pathway ID
##' \item data.pair: the two study names.
##' \item module.type: discordant or concordant modules.
##' }
##' In addition, the elbow plot of -log10(p-value) for each module size will be saved to the filePath.
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step (see the example in function 'merge')
##' dataset.names = c("hb","hs","ht","ha","hi","hl",
##'                   "mb","ms","mt","ma","mi","ml")
##' res_hsa04670 = KEGG_module(mcmc.merge.list, dataset.names,KEGGspecies="hsa",
##'                            KEGGpathwayID="04670",data.pair = c("hs","ms"),
##'                            gene_type = c("discordant"),
##'                            DE_PM_cut = 0.2, minM = 4,maxM = NULL,
##'                            B = 1000, cores = 1,
##'                            search_method = c("Exhaustive"),
##'                            Elbow_plot = T, filePath = getwd())
##' }
KEGG_module = function(mcmc.merge.list,dataset.names,
                       KEGGspecies="hsa",
                       KEGGpathwayID,
                       KEGG.dataGisTopologyG = FALSE,
                       KEGG.dataG2topologyG = NULL,
                       data.pair,
                       gene_type = c("discordant","concordant"),
                       DE_PM_cut = 0.2, minM = 4, maxM = NULL,
                       B = 1000, cores = 1,
                       search_method = c("Exhaustive","SA"),
                       reps_eachM = 1,
                       topG_from_previous=1,
                       Tm0=10,mu=0.95,epsilon=1e-5,N=3000,
                       Elbow_plot = T, filePath = getwd(),
                       seed = 12345, sep = "-"){
  if(minM < 2){
    stop("minM has to be larger than 1.")
  }

  dat1.name = data.pair[[1]]
  dat2.name = data.pair[[2]]
  dat1 = mcmc.merge.list[[match(dat1.name,dataset.names)]]
  dat2 = mcmc.merge.list[[match(dat2.name,dataset.names)]]
  signPM.mat = cbind(apply(dat1,1,mean),apply(dat2,1,mean))

  #match data names and gene names on KEGG topology
  if(KEGG.dataGisTopologyG == TRUE){
    topologyG = rownames(signPM.mat)
  }else if(!is.null(KEGG.dataG2topologyG)){
    topologyG = KEGG.dataG2topologyG[match(rownames(signPM.mat),KEGG.dataG2topologyG[,1]),2]
    na.index = which(is.na(topologyG))
    if(length(na.index != 0)){
      topologyG = topologyG[-na.index]
      signPM.mat = signPM.mat[-na.index,]
    }

  }else if(KEGGspecies == "hsa"){
    map.ls = as.list(org.Hs.eg.db::org.Hs.egALIAS2EG)
    topologyG = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
    na.index = sapply(topologyG, is.null)
    if(sum(na.index) != 0){
      topologyG = topologyG[!na.index]
      signPM.mat = signPM.mat[!na.index,]
    }

  }else if(KEGGspecies == "mmu"){
    map.ls = as.list(org.Mm.eg.db::org.Mm.egALIAS2EG)
    topologyG = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
    na.index = sapply(topologyG, is.null)
    if(sum(na.index) != 0){
      topologyG = topologyG[!na.index]
      signPM.mat = signPM.mat[!na.index,]
    }

  }else if(KEGGspecies == "rno"){
    map.ls = as.list(as.list(org.Rn.eg.db::org.Rn.egALIAS2EG))
    topologyG = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
    na.index = sapply(topologyG, is.null)
    if(sum(na.index) != 0){
      topologyG = topologyG[!na.index]
      signPM.mat = signPM.mat[!na.index,]
    }

  }else if(KEGGspecies == "cel"){
    wormbase = biomaRt::useMart(biomart = "parasite_mart",
                                host = "https://parasite.wormbase.org",
                                port = 443)
    wormbase = useDataset(mart = wormbase, dataset = "wbps_gene")
    map.mat = getBM(attributes = c("entrezgene_name","wormbase_gseq"),
                    filters = "entrezgene_name",
                    values = rownames(signPM.mat),
                    mart = wormbase)
    topologyG = paste0("CELE_",map.mat[match(rownames(signPM.mat),map.mat[,"entrezgene_name"]),"wormbase_gseq"])

    na.index = which(topologyG == "CELE_NA")
    if(length(na.index != 0)){
      topologyG = topologyG[-na.index]
      signPM.mat = signPM.mat[-na.index,]
    }
  }else if(KEGGspecies == "dme"){
    map.ls0 = as.list(org.Dm.eg.db::org.Dm.egALIAS2EG)
    EntrezID = sapply(rownames(signPM.mat),function(g) map.ls0[[g]][[1]])
    map.ls = as.list(org.Dm.eg.db::org.Dm.egFLYBASECG)
    topologyG = paste0("Dmel_",sapply(EntrezID,function(g) ifelse(is.null(g), NA, map.ls[[g]][[1]])))

    na.index = which(topologyG == "Dmel_NA")
    if(length(na.index != 0)){
      topologyG = topologyG[-na.index]
      signPM.mat = signPM.mat[-na.index,]
    }
  }else{
    stop("Please provide mapping between data genes and topology genes when they are of different gene name types and species is not one of 'hsa', 'mmu','rno','cel'or'dme'.")
  }

  rownames(signPM.mat) = topologyG
  adjacent_mat = parseRelation(pathwayID = KEGGpathwayID, keggSpecies = KEGGspecies, sep = sep)
  xmlG = row.names(adjacent_mat)[grep(KEGGspecies,row.names(adjacent_mat))]
  xmlG = gsub(paste0(KEGGspecies,":"),"",xmlG)
  xmlG.ls = lapply(xmlG, function(x){
    strsplit(x,sep)[[1]]
  })

  row.names(adjacent_mat) = colnames(adjacent_mat) = gsub(paste0(KEGGspecies,":"),"",row.names(adjacent_mat))

  mergePMls = lapply(1:length(xmlG.ls), function(x){
    genes = xmlG.ls[[x]]
    cmG = intersect(topologyG,genes)
    if(length(cmG) !=0){
      sub.signPM.mat = matrix(signPM.mat[cmG,],ncol = 2)
      avgPM = apply(sub.signPM.mat, 2, mean)
      return(avgPM)
    }
  })
  names(mergePMls) = xmlG
  mergePMmat = do.call(rbind,mergePMls)

  #discordant/concordant genes definition
  #discordant/concordant genes definition
  if(all(abs(mergePMmat[,1])<=DE_PM_cut | abs(mergePMmat[,2])<=DE_PM_cut)){
    DE_PM_cut = -1
    print(paste0("All genes with DE strength greater than the cutoff value for posterior probability of DE are not connected. Removed the cutoff criterion to consider all ",gene_type," genes regardless of its DE stength."))
  }
  if(gene_type == "discordant"){

    if(sum(mergePMmat[,1]*mergePMmat[,2]<0) == 0){
      stop("No discordant genes are topologically connected. Probably dut to low DE strength in one study. Please check the genePM plot. ")
    }else{
      nodes = unique(rownames(mergePMmat)[which(mergePMmat[,1]*mergePMmat[,2]<0&abs(mergePMmat[,1])>DE_PM_cut&abs(mergePMmat[,2])>DE_PM_cut)])
    }

  }else if(gene_type == "concordant"){

    if(sum(mergePMmat[,1]*mergePMmat[,2]<0) == 0){
      stop("No concordant genes are topologically connected. Probably dut to low DE strength in one study. Please check the genePM plot. ")
    }else{
      nodes = unique(rownames(mergePMmat)[which(mergePMmat[,1]*mergePMmat[,2]>0&abs(mergePMmat[,1])>DE_PM_cut&abs(mergePMmat[,2])>DE_PM_cut)])
    }
  }else{
    stop("gene_type has to be 'discordant' or 'concordant'.")
  }

  undir_adj_mat = adjacent_mat
  for(i in 1:nrow(adjacent_mat)){
    for (j in 1:ncol(adjacent_mat)) {
      undir_adj_mat[i,j] = max(adjacent_mat[i,j],adjacent_mat[j,i])
      undir_adj_mat[j,i] = max(adjacent_mat[i,j],adjacent_mat[j,i])
    }
  }
  g = graph_from_adjacency_matrix(undir_adj_mat,mode="undirected")
  sp.mat <- shortest.paths(g)

  d = degree(g)
  sort.d = sort(d,decreasing = T)[nodes]

  sub.adj.mat = undir_adj_mat[match(nodes,row.names(undir_adj_mat)),
                              match(nodes,colnames(undir_adj_mat))]
  #dim(sub.adj.mat)
  #dim(undir_adj_mat)
  #sum(lower.tri(sub.adj.mat))/length(lower.tri(sub.adj.mat))
  #(sum(lower.tri(undir_adj_mat))-sum(lower.tri(sub.adj.mat)))/(length(lower.tri(undir_adj_mat))-length(lower.tri(sub.adj.mat)))

  if(is.null(maxM)){
    maxM = length(nodes)
  }else{
    maxM = min(length(nodes),maxM)
  }
  module.size = minM:maxM

  search_method = match.arg(search_method)
  if(search_method == "Exhaustive"){
    minG.ls = mclapply(1:length(module.size),function(i){
      m = module.size[i]
      m.combn = combn(x=nodes,m=m)

      ## observed ones:

      asp.m = rep(NA,ncol(m.combn))


      for(j in 1:ncol(m.combn)){
        node.set <- m.combn[,j]
        set.mat <- sp.mat[match(node.set,row.names(sp.mat)),
                          match(node.set,colnames(sp.mat))]
        asp.m[j] <- mean(c(set.mat[lower.tri(set.mat)]))
      }

      (minG = m.combn[,which(asp.m == min(asp.m))])

      set.seed(seed)
      asp.m.perm = rep(NA,B)

      for(b in 1:B){
        permute.set <- sample(xmlG,m)
        permute.mat <- sp.mat[match(permute.set,row.names(sp.mat)),
                              match(permute.set,colnames(sp.mat))]
        asp.m.perm[b] <- mean(c(permute.mat[lower.tri(permute.mat)]))

      }

      p.observed = (sum(asp.m.perm<= min(asp.m)) + 1)/(B+1)
      p.sd = sqrt(p.observed*(1-p.observed)/B)
      null.sp.mean = mean(asp.m.perm)
      null.sp.median = median(asp.m.perm)

      return(list(minG = minG, p.mean = p.observed, p.sd = p.sd,sp = min(asp.m),
                  null.sp.mean = null.sp.mean,
                  null.sp.median = null.sp.median))
    },mc.cores = cores)
  }else{
    minG.ls = list()
    for(i in 1:length(module.size)){
      M = module.size[i]
      print(M)
      if(M>min(module.size)){
        G.ini.list = minG.ls[[i-1]]$top.G
        minG.ls[[i]] = SA_module_M(sp.mat, xmlG, M, nodes, B = B,
                                   G.ini.list=G.ini.list, reps_eachM = reps_eachM,
                                   topG_from_previous=topG_from_previous,
                                   Tm0=Tm0,mu=mu,epsilon=epsilon,
                                   N=N,seed=seed)
      }else{
        minG.ls[[i]] = SA_module_M(sp.mat, xmlG, M, nodes, B = B,
                                   G.ini.list=NULL, reps_eachM = reps_eachM,
                                   topG_from_previous=topG_from_previous,
                                   Tm0=Tm0,mu=mu,epsilon=epsilon,
                                   N=N,seed=seed)
      }
    }
  }
  names(minG.ls) = paste0("minG",module.size)

  #Select best size
  p.mean = sapply(minG.ls, function(x) x[["p.mean"]])
  p.sd = sapply(minG.ls, function(x) x[["p.sd"]])
  #obs.mean.ratio = sapply(minG.ls, function(x) x[["sp"]])/sapply(minG.ls, function(x) x[["null.sp.mean"]])
  obs.median.ratio = sapply(minG.ls, function(x) x[["sp"]])/sapply(minG.ls, function(x) x[["null.sp.median"]])
  obs.median.ratio[is.na(obs.median.ratio)] = 0
  obs.sp = sapply(minG.ls, function(x) x[["sp"]])

  index = which.min(p.mean)
  p.cut = p.mean[index]+2*p.sd[index]
  finalSelect = names(minG.ls)[max(which(p.mean<p.cut & 1:length(p.mean)>=index))]

  KEGGpathwayID_spec = paste0(KEGGspecies,KEGGpathwayID)
  if(Elbow_plot == T){
    names(p.mean) = names(p.sd) = module.size
    pL = sapply(p.mean-2*p.sd, function(x) ifelse(x<0,0,x))
    df = data.frame(logp.observed = -log10(p.mean),
                    logp.max = -log10(pL),
                    logp.min = -log10(p.mean+2*p.sd),
                    size = module.size)
    png(paste0(filePath,"/",KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name,"_",search_method,"_elbow_plot.png"))
    p = ggplot(df, aes(x=size, y=logp.observed)) +
      geom_errorbar(aes(ymin=logp.min, ymax=logp.max), width=.1) +
      geom_line() +
      geom_point() +
      ylim(0,3.1)+
      labs(title = paste0(KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name),y = "-log10(p-value)")
    print(p)
    dev.off()

    df.ratio = data.frame(obs.sp, obs.median.ratio, size = module.size)

    # png(paste0(filePath,"/",KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name,"_",search_method,"_avgSP.png"))
    # p = ggplot(df.ratio, aes(x=size, y=obs.sp)) +
    #   geom_line() +
    #   geom_point()+
    #   labs(title = paste0(KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name),y = "module average shortest path value")
    # print(p)
    # dev.off()

    # png(paste0(filePath,"/",KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name,"_",search_method,"avgSP_ratio_obs_to_median.png"))
    # p = ggplot(df.ratio, aes(x=size, y=obs.median.ratio)) +
    #   geom_line() +
    #   geom_point() +
    #   labs(title = paste0(KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name),y = "average module shortest path/median(null average shortest path)")
    # print(p)
    # dev.off()

  }
  return(list(minG.ls=minG.ls,bestSize = finalSelect,
              mergePMmat = mergePMmat,
              KEGGspecies = KEGGspecies,
              KEGGpathwayID = KEGGpathwayID,
              data.pair = data.pair,
              module.type = gene_type))
}


