##' The \code{KEGG_module_topology_plot} is a function to highlght module location on KEGG pathway
##' topology.
##' @title KEGG topology plot with highlighted module nodes.
##' @param res_KEGG_module: a result list from function \code{KEGG_module}.
##' @param which_to_draw: either a numeric vector indicating module sizes of interest to highlight or "all"
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
  module.size = as.numeric(gsub("minG","",names(minG.ls)))
  module.type = res_KEGG_module$module.type

  mergePMmat = res_KEGG_module$mergePMmat
  KEGGspecies = res_KEGG_module$KEGGspecies
  KEGGpathwayID = res_KEGG_module$KEGGpathwayID
  KEGGpathwayID_spec = paste0(KEGGspecies,KEGGpathwayID)
  data.pair = res_KEGG_module$data.pair
  dat1.name = data.pair[1]
  dat2.name = data.pair[2]

  if(which_to_draw[[1]] == "all"){
    which_to_draw_index = 1:length(minG.ls)
  }else if(!is.numeric(which_to_draw)){
    stop("which_to_draw should be 'all' or a numeric vector")
  }else{
    which_to_draw_index = match(which_to_draw,module.size)
  }
  orig.path = getwd()
  setwd(filePath)
  for (j in 1:length(which_to_draw_index)) {
    index = which_to_draw_index[j]
    topologyG = minG.ls[[index]]$minG
    if(is.null(dim(topologyG))){
      signPM.mat = mergePMmat[topologyG,]
      row.names(signPM.mat) = sapply(row.names(signPM.mat), function(x) strsplit(x,"_")[[1]][1])

      res = pathview(gene.data = signPM.mat, pathway.id = KEGGpathwayID,
                     species = KEGGspecies, out.suffix = "", kegg.native = T,
                     key.pos = "bottomright", map.null=T,cex = 0.15)

      file.rename(paste(KEGGpathwayID_spec,"..multi.png",sep=""),
                  paste(KEGGpathwayID_spec,"_",dat1.name,"_",dat2.name,"_",names(minG.ls)[index],"_",module.type,".png",sep=""))
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
                    paste(KEGGpathwayID_spec,"_",dat1.name,"_",dat2.name,"_",names(minG.ls)[index],"_",i,"_",module.type,".png",sep=""))
        file.remove(paste(KEGGpathwayID_spec,".xml",sep=""))
        file.remove(paste(KEGGpathwayID_spec,".png",sep=""))

      }
    }
  }
  setwd(orig.path)
}
