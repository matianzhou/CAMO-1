##' Read raw data: microarray & RNAseq 
##' The \code{read.rawData} is function to read raw data
##' @title Read raw data: microarray or RNAseq.  
##' @param data.file: the file name of the data to be read.
##' @param sep: the field separator character. Values on each line of
##' the file are separated by this character. If \code{sep = ""}, the
##' separator is 'white space', that is one or more spaces, tabs,
##' newlines or carriage returns. The default is \code{,}. See
##' \link{read.csv}.
##' @param quote: the set of quoting characters. To disable quoting
##' altogether, use \code{quote = ""}.  See \link{scan} for the
##' behaviour on quotes embedded in quotes. Quoting is only
##' considered for columns read as character, which is all of them
##' unless colClasses is specified.  The default is \code{\"}. See
##' \link{read.csv} 
##' @param header: a logical value indicating whether the file contains
##' the names of the variables as its first line. If missing, the
##' value is determined from the file format: \code{header} is set to
##' \code{TRUE} if and only if the first row contains one fewer field
##' than the number of columns. The default is \code{TRUE}. See
##' \link{read.csv}.

##' @return the data matrix.
##' @export
##' @examples
##' \dontrun{
##' data(hb)
##' write.csv(data,file="rawData.csv")
##' data <- read.rawData(data.file="rawData.csv")
##' }

read.rawData <- function(data.file, sep=",", quote='"', header=T) {
  #input either csv or txt file, gene on rows and sample on columns
    data <- as.matrix(read.csv(data.file, sep=",", quote=quote, 
                               header=header, row.names=1))
    check.rawData(data)  # check the expression data 
    return(data)
}

##' Read group label data
##' The \code{read.groupData} is function to read group label data
##' @title Read group label data.
##' @param group.file: the file name of the group data to be read.
##' @param sep: the field separator character. Values on each line of
##' the file are separated by this character. If \code{sep = ""}, the
##' separator is 'white space', that is one or more spaces, tabs,
##' newlines or carriage returns. The default is \code{,}. See
##' \link{read.csv}.
##' @param quote: the set of quoting characters. To disable quoting
##' altogether, use \code{quote = ""}.  See \link{scan} for the
##' behaviour on quotes embedded in quotes. Quoting is only
##' considered for columns read as character, which is all of them
##' unless colClasses is specified.  The default is \code{\"}. See
##' \link{read.csv} 
##' @param header: a logical value indicating whether the file contains
##' the names of the variables as its first line. If missing, the
##' value is determined from the file format: \code{header} is set to
##' \code{TRUE} if and only if the first row contains one fewer field
##' than the number of columns. The default is \code{TRUE}. See
##' \link{read.csv}.

##' @return the group label vector.
##' @export
##' @examples
##' \dontrun{
##' data(hb)
##' write.csv(group,file="group.csv")
##' group <- read.groupData(group.file = "group.csv")
##' }

read.groupData <- function(group.file, sep=",", quote='"', header=T) {
  #input either csv or txt file, samples on rows, one column containing the class label
  group <- as.factor(read.csv(group.file, sep=",", quote=quote, 
                              header=header, row.names=1)[,1])
  check.groupData(group)  # check the group data
  return(group)
}
  
##' Read p-value data
##' The \code{read.pData} is function to read p-value data
##' @title Read p-value data.
##' @param p.file: the file name of the group data to be read.
##' @param sep: the field separator character. Values on each line of
##' the file are separated by this character. If \code{sep = ""}, the
##' separator is 'white space', that is one or more spaces, tabs,
##' newlines or carriage returns. The default is \code{,}. See
##' \link{read.csv}.
##' @param quote: the set of quoting characters. To disable quoting
##' altogether, use \code{quote = ""}.  See \link{scan} for the
##' behaviour on quotes embedded in quotes. Quoting is only
##' considered for columns read as character, which is all of them
##' unless colClasses is specified.  The default is \code{\"}. See
##' \link{read.csv} 
##' @param header: a logical value indicating whether the file contains
##' the names of the variables as its first line. If missing, the
##' value is determined from the file format: \code{header} is set to
##' \code{TRUE} if and only if the first row contains one fewer field
##' than the number of columns. The default is \code{TRUE}. See
##' \link{read.csv}.

##' @return the group label vector.
##' @export
##' @examples
##' \dontrun{
##' data(hb)
##' summaryDE <- indDE(data=data,group=group,data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' pData <- summaryDE[,c(3,1)]
##' write.csv(pData,file="pData.csv")
##' pData <- read.pData(p.file="pData.csv")
##' }

read.pData <-  function(p.file, sep=",",quote='"', header=T){
  #input either csv or txt file, gene on rows, two columns: 1st: 2-sided p-value, 2nd: logFC (or equivalent effect size)
  pData <- as.matrix(read.csv(p.file, sep=",", quote=quote, 
                              header=header, row.names=1))
  check.pData(pData)  # check the pvalue data 
  return(pData)
} 

