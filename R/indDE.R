##' Differential expression analysis for individual data
##' The \code{indDE} is function to perform differential expression
##' analysis for individual data
##' @title Differential expression analysis for individual data.
##' @param data: the raw expression data.
##' @param group: the group label.
##' @param data.type: either "microarray" or "RNAseq"
##' @param case.label: label for the case group.
##' @param control.label: label for the control group.

##' @return summary consisting of log fold change, lfc standard error,
##' p-value, q-value.
##' @export
##' @examples
##' \dontrun{
##' data(hb)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                    case.label="2", ctrl.label="1")
##' }
indDE <- function(data, group, data.type, case.label, ctrl.label){
  ## data is a numeric matrix, group is a factor

  check.compatibility(data, group, case.label, ctrl.label)

  if(data.type=="microarray") {
    #implement limma
    group <- relevel(group,ref=ctrl.label)
    design <-model.matrix(~group)
    fit <-lmFit(data, design)
    ebFit<-eBayes(fit)
    out.table <- topTable(ebFit,coef=2, number=Inf, sort.by='none')
    log2FC <- out.table$logFC
    lfcSE <- sqrt(ebFit$s2.post) * fit$stdev.unscaled[,2]
    p <- as.numeric(out.table$P.Value)
    q <- p.adjust(p,method="BH")
    res <- data.frame(logFC=log2FC,lfcSE = lfcSE,
                          pvalue = p, qvalue= q)
    rownames(res) <- rownames(data)
  }else if(data.type=="RNAseq") {
    #implement DESeq2
    group <- relevel(group,ref=ctrl.label)
    design <-model.matrix(~ group)  # design matrix
    colData <- data.frame(group=group)
    colnames(colData) <- colnames(design)[-1]
    ddsMat <- DESeqDataSetFromMatrix(countData = data,
                                     colData = colData,
                                     design = as.formula(
                                       paste(" ~ ",paste(colnames(colData), collapse=" + ")
                                       ) )  )
    ddsMat <- DESeq(ddsMat)
    res <- results(ddsMat,contrast=c(colnames(colData)[1],levels(group)[2],
                                     levels(group)[1]))
    log2FC <- as.numeric(res$log2FoldChange)
    lfcSE <- as.numeric(res$lfcSE)
    p <- as.numeric(res$pvalue)
    q <- p.adjust(p,method="BH")
    res <- data.frame(logFC=log2FC,lfcSE = lfcSE,
                          pvalue = p, qvalue= q)
    rownames(res) <- rownames(data)
  }else{
    stop("data.type has to be 'microarray' or 'RNAseq'.")
  }
  return(res)

}
