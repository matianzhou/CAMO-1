##' KEGG pathway topology with highlighted module nodes
##' The \code{KEGG_module_topology_plot} is a function to highlght module location on KEGG pathway
##' topology.
##' @title KEGG pathway topology with highlighted module nodes.
##' @param res_KEGG_module: a result list from function \code{KEGG_module}.
##' @param which_to_draw: either a numeric vector indicating which modules to highlight or "all"
##' which will show all module size scenarios in the res_KEGG_module provided.
##' @param filePath: the path to save the elbow plot. Default is the current working directory.
##' @return KEGG pathway topology plots with different module highlighted will be saved as .png
##' files in the filePath provided.
##' @export
##' @examples
##' \dontrun{
##' #res_KEGG_module from the KEGG_module step (see the example in function 'KEGG_module')
##' res = KEGG_module_topology_plot(res_KEGG_module,which_to_draw = c(4,8,9))
##' }
KEGG_module_topology_plot = function(res_KEGG_module,which_to_draw = "all",filePath = getwd()){
  minG.ls = res_KEGG_module$minG.ls
  mergePMmat = res_KEGG_module$mergePMmat
  KEGGspecies = res_KEGG_module$KEGGspecies
  KEGGpathwayID = res_KEGG_module$KEGGpathwayID
  KEGGpathwayID_spec = paste0(KEGGspecies,KEGGpathwayID)
  data.pair = res_KEGG_module$data.pair
  dat1.name = data.pair[1]
  dat2.name = data.pair[2]

  if(which_to_draw[[1]] == "all"){
    which_to_draw = 1:length(minG.ls)
  }else if(!is.numeric(which_to_draw)){
    stop("which_to_draw should be 'all' or a numeric vector")
  }
  orig.path = getwd()
  setwd(filePath)
  for (j in 1:length(which_to_draw)) {
    index = which_to_draw[j]
    topologyG = minG.ls[[index]]$minG
    if(is.null(dim(topologyG))){
      signPM.mat = mergePMmat[topologyG,]
      row.names(signPM.mat) = sapply(row.names(signPM.mat), function(x) strsplit(x,"_")[[1]][1])

      res = pathview(gene.data = signPM.mat, pathway.id = KEGGpathwayID,
                     species = KEGGspecies, out.suffix = "", kegg.native = T,
                     key.pos = "bottomright", map.null=T,cex = 0.15)

      file.rename(paste(KEGGpathwayID_spec,"..multi.png",sep=""),
                  paste(KEGGpathwayID_spec,"_",dat1.name,"_",dat2.name,"_",names(minG.ls)[index],".png",sep=""))
      file.remove(paste(KEGGpathwayID_spec,".xml",sep=""))
      file.remove(paste(KEGGpathwayID_spec,".png",sep=""))


    }else{
      for (i in 1:ncol(topologyG)) {
        topologyG0 = topologyG[,i]
        signPM.mat = mergePMmat[topologyG0,]
        row.names(signPM.mat) = sapply(row.names(signPM.mat), function(x) strsplit(x,"_")[[1]][1])

        res = pathview(gene.data = signPM.mat, pathway.id = KEGGpathwayID,
                       species = KEGGspecies, out.suffix = "", kegg.native = T,
                       key.pos = "bottomright", map.null=T,cex = 0.15)
        file.rename(paste(KEGGpathwayID_spec,"..multi.png",sep=""),
                    paste(KEGGpathwayID_spec,"_",dat1.name,"_",dat2.name,"_",names(minG.ls)[index],"_",i,".png",sep=""))
        file.remove(paste(KEGGpathwayID_spec,".xml",sep=""))
        file.remove(paste(KEGGpathwayID_spec,".png",sep=""))

      }
    }
  }
  setwd(orig.path)
}
