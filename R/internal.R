##########################
###### Check functions ###
##########################

check.rawData <- function(data){
  G <- nrow(data)
  N <- ncol(data)
  if(!is.numeric(data)) {
    stop("expression data not numeric")
  }
  if(is.null(row.names(data))) {
    stop("gene symbol missing")
  }
  if(N <= 3) {
    stop("too few samples")
  }
  if(sum(duplicated(row.names(data)))>0){
  	stop("duplicate gene symbols")
  }
}

check.groupData <- function(group){
  l <- nlevels(group)
  if( l != 2 ) {
    stop("not a two-class comparison")
  }
}

check.pData <- function(pData){
  G <- nrow(pData)
  if(!is.numeric(pData)) {
    stop("pvalue data not numeric")
  }
  if(is.null(row.names(pData))) {
    stop("gene symbol missing")
  }
  if(ncol(pData) <2) {
    stop("missing either p-value or effect size")
  }
}


check.compatibility <- function(data, group, case.label, ctrl.label){
  G <- nrow(data)
  N1 <- ncol(data)
  N2 <- length(group)
  if(N1 != N2) {
    stop("expression data and class label have unmatched sample size")
  }
  if(!all(group %in% c(case.label,ctrl.label))){
    stop("including class labels other than the case and control")
  }
  if(sum(group==case.label) <= 1 ||sum(group==ctrl.label) <= 1){
    stop("not enough samples in either case or control group")
  }
}

##########################
###### BayesP part ###
##########################

PtoZ <- function(p2, lfc) {
  sgn <- sign(lfc)
  z <- ifelse(sgn>0, qnorm(p2/2), qnorm(p2/2,lower.tail = F))
  return(z)
}

SelectGamma <- function(p){
  ## Gamma is DE proportion  = 1-pi0
  m <- length(p)
  lambda <- seq(0,0.95,by=0.01)
  pi0 <- sapply(lambda, function(x) sum(p>x)/(m*(1-x))  )

  # fit a natural cubic spline
  library(splines)
  dat <- data.frame(pi0=pi0, lambda=lambda)
  lfit <- lm(pi0 ~ ns(lambda, df = 3), data=dat)
  pi0hat <- predict(lfit, data.frame(lambda=1))
  gamma <- 1 - pi0hat
  return(gamma)
}

##################################
###### ACS scores ################
##################################

CSY <- function(d1,d2,index1){
  Sens <- calcSensC(d1[index1,],d2[index1,])
  Spec <- calcSpecC(d1[-index1,],d2[-index1,])

  CS <- Sens + Spec - 1
  if(is.nan(CS)) {
    CS <- 0
  }
  return(CS)
}

CSF <- function(d1,d2,index1,index2){
  Sens <- calcSensC(d1[index1,],d2[index1,])
  Prec <- calcPrecC(d1[index2,],d2[index2,])

  CS <- (2*Sens*Prec)/(Sens+Prec)
  if(is.nan(CS)) {
    CS <- 0
  }
  return(CS)
}

CSG <- function(d1,d2,index1,index2){
  Sens <- calcSensC(d1[index1,],d2[index1,])
  Spec <- calcSpecC(d1[-index1,],d2[-index1,])

  CS <- sqrt(Sens*Spec)
  if(is.nan(CS)) {
    CS <- 0
  }
  return(CS)
}


ECSY <- function(d1,d2){
  ESens <- calcESensC(d1,d2)
  ESpec <- calcESpecC(d1,d2)

  ECS <- ESens + ESpec - 1
  if(is.nan(ECS)) {
    ECS <- 0
  }
  return(ECS)
}

ECSF <- function(d1,d2){
  ESens <- calcESensC(d1,d2)
  EPrec <- calcEPrecC(d1,d2)

  ECS <- (2*ESens*EPrec)/(ESens+EPrec)
  if(is.nan(ECS)) {
    ECS <- 0
  }
  return(ECS)
}

ECSG <- function(d1,d2){
  ESens <- calcESensC(d1,d2)
  ESpec <- calcESpecC(d1,d2)

  ECS <- sqrt(ESens*ESpec)
  if(is.nan(ECS)) {
    ECS <- 0
  }
  return(ECS)
}

permCSY <- function(d1,d2){
  Sens <- calcSensC(d1,d2)
  Spec <- calcSpecC(d1,d2)

  permCS <- Sens + Spec - 1
  if(is.nan(permCS)) {
    permCS <- 0
  }
  return(permCS)
}

permCSF <- function(d1,d2){
  Sens <- calcSensC(d1,d2)
  Prec <- calcPrecC(d1,d2)

  permCS <- (2*Sens*Prec)/(Sens+Prec)
  if(is.nan(permCS)) {
    permCS <- 0
  }
  return(permCS)
}

permCSG <- function(d1,d2){
  Sens <- calcSensC(d1,d2)
  Spec <- calcSpecC(d1,d2)

  permCS <- sqrt(Sens*Spec)
  if(is.nan(permCS)) {
    permCS <- 0
  }
  return(permCS)
}


CS <- function(dat1,dat2,deIndex1,deIndex2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    CS <- CSY(dat1,dat2,deIndex1)
  } else if(measure=="Fmeasure"){
    CS <- CSF(dat1,dat2,deIndex1,deIndex2)
  } else if(measure=="geo.mean"){
    CS <- CSG(dat1,dat2,deIndex1)
  }
  return(CS)
}


ECS <- function(dat1,dat2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    ECS <- ECSY(dat1,dat2)
  } else if(measure=="Fmeasure"){
    ECS <- ECSF(dat1,dat2)
  } else if(measure=="geo.mean"){
    ECS <- ECSG(dat1,dat2)
  }
  return(ECS)
}

permCS <- function(dat1,dat2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    pemrCS <- permCSY(dat1,dat2)
  } else if(measure=="Fmeasure"){
    permCS <- permCSF(dat1,dat2)
  } else if(measure=="geo.mean"){
    permCS <- permCSG(dat1,dat2)
  }
  return(permCS)
}


##################################
###### ADS scores ################
##################################


DSY <- function(d1,d2,index1){
  Sens <- calcSensD(d1[index1,],d2[index1,])
  Spec <- calcSpecD(d1[-index1,],d2[-index1,])

  DS <- Sens + Spec - 1
  if(is.nan(DS)) {
    DS <- 0
  }
  return(DS)
}

DSF <- function(d1,d2,index1,index2){
  Sens <- calcSensD(d1[index1,],d2[index1,])
  Prec <- calcPrecD(d1[index2,],d2[index2,])

  DS <- (2*Sens*Prec)/(Sens+Prec)
  if(is.nan(DS)) {
    DS <- 0
  }
  return(DS)
}

DSG <- function(d1,d2,index1,index2){
  Sens <- calcSensD(d1[index1,],d2[index1,])
  Spec <- calcSpecD(d1[-index1,],d2[-index1,])

  DS <- sqrt(Sens*Spec)
  if(is.nan(DS)) {
    DS <- 0
  }
  return(DS)
}


EDSY <- function(d1,d2){
  ESens <- calcESensD(d1,d2)
  ESpec <- calcESpecD(d1,d2)

  EDS <- ESens + ESpec - 1
  if(is.nan(EDS)) {
    EDS <- 0
  }
  return(EDS)
}

EDSF <- function(d1,d2){
  ESens <- calcESensD(d1,d2)
  EPrec <- calcEPrecD(d1,d2)

  EDS <- (2*ESens*EPrec)/(ESens+EPrec)
  if(is.nan(EDS)) {
    EDS <- 0
  }
  return(EDS)
}

EDSG <- function(d1,d2){
  ESens <- calcESensD(d1,d2)
  ESpec <- calcESpecD(d1,d2)

  EDS <- sqrt(ESens*ESpec)
  if(is.nan(EDS)) {
    EDS <- 0
  }
  return(EDS)
}

permDSY <- function(d1,d2){
  Sens <- calcSensD(d1,d2)
  Spec <- calcSpecD(d1,d2)

  permDS <- Sens + Spec - 1
  if(is.nan(permDS)) {
    permDS <- 0
  }
  return(permDS)
}

permDSF <- function(d1,d2){
  Sens <- calcSensD(d1,d2)
  Prec <- calcPrecD(d1,d2)

  permDS <- (2*Sens*Prec)/(Sens+Prec)
  if(is.nan(permDS)) {
    permDS <- 0
  }
  return(permDS)
}

permDSG <- function(d1,d2){
  Sens <- calcSensD(d1,d2)
  Spec <- calcSpecD(d1,d2)

  permDS <- sqrt(Sens*Spec)
  if(is.nan(permDS)) {
    permDS <- 0
  }
  return(permDS)
}


DS <- function(dat1,dat2,deIndex1,deIndex2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    DS <- DSY(dat1,dat2,deIndex1)
  } else if(measure=="Fmeasure"){
    DS <- DSF(dat1,dat2,deIndex1,deIndex2)
  } else if(measure=="geo.mean"){
    DS <- DSG(dat1,dat2,deIndex1)
  }
  return(DS)
}


EDS <- function(dat1,dat2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    EDS <- EDSY(dat1,dat2)
  } else if(measure=="Fmeasure"){
    EDS <- EDSF(dat1,dat2)
  } else if(measure=="geo.mean"){
    EDS <- EDSG(dat1,dat2)
  }
  return(EDS)
}

permDS <- function(dat1,dat2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    pemrDS <- permDSY(dat1,dat2)
  } else if(measure=="Fmeasure"){
    permDS <- permDSF(dat1,dat2)
  } else if(measure=="geo.mean"){
    permDS <- permDSG(dat1,dat2)
  }
  return(permDS)
}



##################################
#### Global (expected value from marginal,
#### permutate genes to get p-value)
##################################


perm_global <- function(dat1,dat2,measure="Fmeasure",B){
  G <- nrow(dat1)

  #rawCS <- CS(dat1,dat2,deIndex1,deIndex2,measure)
  #expCS <- ECS(dat1,dat2,measure)
  #rawDS <- DS(dat1,dat2,deIndex1,deIndex2,measure)
  #expDS <- EDS(dat1,dat2,measure)

  out <- matrix(NA,B,4)
  colnames(out) <- c("permCS","permECS","permDS","permEDS")

  for(b in 1:B){
    dat1perm <- dat1[sample(1:G,G,replace = F),]
    dat2perm <- dat2[sample(1:G,G,replace = F),]
    out[b,"permCS"] <- permCS(dat1perm,dat2perm,measure)
    out[b,"permECS"] <- ECS(dat1perm,dat2perm,measure)
    out[b,"permDS"] <- permDS(dat1perm,dat2perm,measure)
    out[b,"permEDS"] <- EDS(dat1perm,dat2perm,measure)
  }
  return(out)
}

ACS_global <- function(dat1,dat2,deIndex1,deIndex2,
                       measure="Fmeasure"){
  cs <- CS(dat1,dat2,deIndex1,deIndex2,measure)
  ecs <- ECS(dat1,dat2,measure)
  acs <- (cs - ecs)/(1-ecs)
  return(acs)
}

pACS_global <- function(dat1,dat2,deIndex1,deIndex2,
                        measure="Fmeasure",acs,permOut){

  permcs <- permOut[,"permCS"]
  permecs <- permOut[,"permECS"]
  permacs <- (permcs - permecs)/(1-permecs)

  p_acs <- (sum(permacs>=acs) + 1)/(length(permacs)+1)
  return(p_acs)
}

ADS_global <- function(dat1,dat2,deIndex1,deIndex2,
                       measure="Fmeasure"){
  ds <- DS(dat1,dat2,deIndex1,deIndex2,measure)
  eds <- EDS(dat1,dat2,measure)
  ads <- (ds - eds)/(1-eds)
  return(ads)
}

pADS_global <- function(dat1,dat2,deIndex1,deIndex2,
                        measure="Fmeasure",ads,permOut){
  permds <- permOut[,"permDS"]
  permeds <- permOut[,"permEDS"]
  permads <- (permds - permeds)/(1-permeds)

  p_ads <- (sum(permads>=ads) + 1)/(length(permads)+1)
  return(p_ads)
}


##################################
#### Pathway (expected value from global,
#### permute genes to get p-value)
##################################


margin_pathway <-  function(dat1,dat2,
                            select.pathway.list,
                            measure="Fmeasure"){
  select.pathways <- names(select.pathway.list)
  data_genes <- rownames(dat1)
  pathway.size <- sapply(select.pathway.list,function(x) {
    length(intersect(data_genes,x))})
  K <- length(select.pathways)
  G <- nrow(dat1)

  out <- matrix(NA,nrow=K,ncol=2)
  rownames(out) <- select.pathways
  colnames(out) <- c("ECS","EDS")

  R <- 20 ##fairly enough

  for(k in 1:K){
    #print(k)
    ecsk <- edsk <- rep(NA,R)
    pathsizek <- pathway.size[k]
    for(j in 1:R){
      index <- sample(1:G,pathsizek,replace=F)
      dat1.select <- dat1[index,]
      dat2.select <- dat2[index,]
      ecsk[j] <- ECS(dat1.select,dat2.select,measure)
      edsk[j] <- EDS(dat1.select,dat2.select,measure)
    }
    out[k,"ECS"] <- mean(ecsk)
    out[k,"EDS"] <- mean(edsk)
  }

  return(out)

}

perm_pathway <- function(dat1,dat2,
                         select.pathway.list,
                         measure="Fmeasure",B,parallel=F,n.cores=4){

  select.pathways <- names(select.pathway.list)
  data_genes <- rownames(dat1)
  pathway.size <- sapply(select.pathway.list,function(x) {
    length(intersect(data_genes,x))})
  K <- length(select.pathways)
  G <- nrow(dat1)

  #out <- array(1,dim=c(B,K,4),dimnames=
  #list(1:B,select.pathways,
  #c("permCS","permECS","permDS","permEDS")))

  #rawCS <- CS(dat1,dat2,deIndex1,deIndex2,measure)
  #expCS <- ECS(dat1,dat2,measure)
  #rawDS <- DS(dat1,dat2,deIndex1,deIndex2,measure)
  #expDS <- EDS(dat1,dat2,measure)

  out <- array(1,dim=c(B,K,2),dimnames=
                 list(1:B,select.pathways,c("permCS","permDS")))

  for(k in 1:K){
    #print(k)
    pathsizek <- pathway.size[k]
    if(parallel == T){
      permFunc = function(b){
        dat1perm <- dat1[sample(1:G,G,replace = F),]
        dat2perm <- dat2[sample(1:G,G,replace = F),]
        index <- sample(1:G,pathsizek,replace=F)
        dat1perm.select <- dat1perm[index,]
        dat2perm.select <- dat2perm[index,]
        permCS_res <- permCS(dat1perm.select,dat2perm.select,measure)
        permDS_res <- permDS(dat1perm.select,dat2perm.select,measure)
        return(list(permCS_res = permCS_res, permDS_res = permDS_res))
      }
      out.ls = mclapply(1:B, permFunc, mc.cores = n.cores)
      for(b in 1:B){
        out[b,k,"permCS"] <- out.ls[[b]]$permCS_res
        out[b,k,"permDS"] <- out.ls[[b]]$permDS_res
      }
    }else{
      for(b in 1:B){
        dat1perm <- dat1[sample(1:G,G,replace = F),]
        dat2perm <- dat2[sample(1:G,G,replace = F),]
        index <- sample(1:G,pathsizek,replace=F)
        dat1perm.select <- dat1perm[index,]
        dat2perm.select <- dat2perm[index,]

        out[b,k,"permCS"] <- permCS(dat1perm.select,dat2perm.select,measure)
        #out[b,k,"permECS"] <- ECS(dat1perm.select,dat2perm.select,measure)
        out[b,k,"permDS"] <- permDS(dat1perm.select,dat2perm.select,measure)
        #out[b,k,"permEDS"] <- EDS(dat1perm.select,dat2perm.select,measure)

      }
    }
  }
  return(out)
}

ACS_pathway <- function(dat1,dat2,deIndex1,deIndex2,
                        select.pathway.list,
                        measure="Fmeasure",marginOut){

  select.pathways <- names(select.pathway.list)
  data_genes <- rownames(dat1)
  pathway.size <- sapply(select.pathway.list,function(x) {
    length(intersect(data_genes,x))})
  K <- length(select.pathways)
  G <- nrow(dat1)

  acs <- rep(NA,K)
  names(acs) <- select.pathways

  for(k in 1:K){
    path_genek <- select.pathway.list[[k]]
    genek <- intersect(path_genek,data_genes)
    dat1_k <- dat1[genek,]
    dat2_k <- dat2[genek,]

    if(length(intersect(names(deIndex1), genek))<=3 ){
      deIndex1_k <- 1:nrow(dat1_k)
    } else {
      deIndex1_k <- which(rownames(dat1_k)%in%intersect(names(deIndex1), genek))
    }

    if(length(intersect(names(deIndex2), genek))<=3 ){
      deIndex2_k <- 1:nrow(dat2_k)
    } else {
      deIndex2_k <- which(rownames(dat2_k)%in%intersect(names(deIndex2), genek))
    }

    cs <- CS(dat1_k,dat2_k,deIndex1_k,deIndex2_k,measure)
    ecs <- marginOut[k,"ECS"]
    acs[k] <- (cs - ecs)/(1-ecs)
  }
  return(acs)
}

pACS_pathway <- function(dat1,dat2,deIndex1,deIndex2,
                         select.pathway.list,
                         measure="Fmeasure",acs,permOut,marginOut){

  select.pathways <- names(select.pathway.list)
  K <- length(select.pathways)

  p_acs <- rep(NA,K)
  names(p_acs) <- select.pathways

  for(k in 1:K){

    permcs <- permOut[,k,"permCS"]
    ecs <- marginOut[k,"ECS"]
    #permecs <- permOut[,k,"permECS"]
    permacs <- (permcs - ecs)/(1-ecs)

    p_acs[k] <- (sum(permacs>=acs[k]) + 1)/(length(permacs)+1)
  }

  return(p_acs)

}


ADS_pathway <- function(dat1,dat2,deIndex1,deIndex2,
                        select.pathway.list,
                        measure="Fmeasure",marginOut){

  select.pathways <- names(select.pathway.list)
  data_genes <- rownames(dat1)
  pathway.size <- sapply(select.pathway.list,function(x) {
    length(intersect(data_genes,x))})
  K <- length(select.pathways)
  G <- nrow(dat1)

  ads <- rep(NA,K)
  names(ads) <- select.pathways

  for(k in 1:K){
    path_genek <- select.pathway.list[[k]]
    genek <- intersect(path_genek,data_genes)
    dat1_k <- dat1[genek,]
    dat2_k <- dat2[genek,]

    if(length(intersect(names(deIndex1), genek))<=3 ){
      deIndex1_k <- 1:nrow(dat1_k)
    } else {
      deIndex1_k <- which(rownames(dat1_k)%in%intersect(names(deIndex1), genek))
    }

    if(length(intersect(names(deIndex2), genek))<=3 ){
      deIndex2_k <- 1:nrow(dat2_k)
    } else {
      deIndex2_k <- which(rownames(dat2_k)%in%intersect(names(deIndex2), genek))
    }

    ds <- DS(dat1_k,dat2_k,deIndex1_k,deIndex2_k,measure)
    eds <- marginOut[k,"EDS"]
    ads[k] <- (ds - eds)/(1-eds)
  }
  return(ads)
}

pADS_pathway <- function(dat1,dat2,deIndex1,deIndex2,
                         select.pathway.list,
                         measure="Fmeasure",ads,permOut,marginOut){

  select.pathways <- names(select.pathway.list)
  K <- length(select.pathways)

  p_ads <- rep(NA,K)
  names(p_ads) <- select.pathways

  for(k in 1:K){

    permds <- permOut[,k,"permDS"]
    #permeds <- permOut[,k,"permEDS"]
    eds <- marginOut[k,"EDS"]
    permads <- (permds - eds)/(1-eds)

    p_ads[k] <- (sum(permads>=ads[k]) + 1)/(length(permads)+1)

  }

  return(p_ads)

}

##########################
##Pathway enrich analysis#
##########################

gsa.fisher <- function(x, background, pathway) {
  ####x is the list of query genes
  ####backgroud is a list of background genes that query genes from
  ####pathway is a list of different pathway genes
  count_table<-matrix(0,2,2)
  x<-toupper(x)
  background<-toupper(background)
  index<-which(toupper(background) %in% toupper(x)==FALSE)
  background_non_gene_list<-background[index]
  x<-toupper(x)
  pathway<-lapply(pathway,function(x) intersect(toupper(background),toupper(x)))
  get.fisher <- function(path) {
    res <- NA
    ####in the gene list and in the pathway
    count_table[1,1]<-sum(x %in% path)
    #count_table[1,1]<-sum(is.na(charmatch(x,path))==0)
    ####in the gene list but not in the pathway
    count_table[1,2]<-length(x)-count_table[1,1]
    ####not in the gene list but in the pathway
    count_table[2,1]<-sum(background_non_gene_list%in% path)
    ####not in the gene list and not in the pathway
    count_table[2,2]<-length(background_non_gene_list)-count_table[2,1]
    matched_gene<-x[x %in% path]
    match_num<-length(matched_gene)
    overlap_info<-array(0,dim=4)
    names(overlap_info)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
    overlap_info[1]=count_table[1,1]
    overlap_info[2]=count_table[1,2]
    overlap_info[3]=count_table[2,1]
    overlap_info[4]=count_table[2,2]
    if(length(count_table)==4){
      res <- fisher.test(count_table, alternative="greater")$p}
    return(list(p_value=res,match_gene=matched_gene,match_num=match_num,
                fisher_table=overlap_info))
  }
  p_val<-rep(0,length(pathway))

  match_gene_list<-list(length(pathway))
  match_gene <- array(0,dim=length(pathway))
  num1<-array(0,dim=length(pathway))
  num2<-matrix(0,nrow=length(pathway),ncol=4)
  colnames(num2)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
  for(i in 1:length(pathway)){
    result<-get.fisher(pathway[[i]])
    p_val[i]<-result$p_value
    match_gene_list[[i]]<-result$match_gene
    match_gene[i]<-paste(match_gene_list[[i]],collapse="/")
    num1[i]<-result$match_num
    num2[i,]<-result$fisher_table
  }
  names(p_val) <- names(pathway)
  q_val <- p.adjust(p_val, "BH")

  summary<-data.frame(pvalue=p_val,
                      qvalue=q_val,
                      DE_in_Set=num2[,1],
                      DE_not_in_Set=num2[,2],
                      NonDE_in_Set=num2[,3],
                      NonDE_not_in_Set=num2[,4])

  a<-format(summary,digits=3)

  return(a)
}

fisher <- function(x){
  n <- length(x)
  y <- -2*log(x)
  Tf <- sum(y)
  return(1-pchisq(Tf,2*n))
}


##########################
###### SA functions ######
##########################

## energy function

E_tot <- function(delta.mat,a,delta.est){
  ## vector "a" of length n
  ## vector "delta_est" of length K+1: start from theta_0, then ordered from k=1 to K
  n <- nrow(delta.mat)
  K <- length(unique(a))
  #theta_0 <- delta.est[1]
  theta_0 <- 0
  E <- sum(sapply(1:n, function(x) {
    sum(sapply(1:n, function(y){
      if(a[x]==a[y]){
        (delta.mat[x,y] - delta.est[as.character(a[x])])^2
      } else{
        (delta.mat[x,y] - theta_0)^2
      }
    },simplify=T))
  }, simplify=T))

  return(E)
}

## Estimate of delta.est (the means)

Est_mean <- function(delta.mat,a){
  # the first element is always the off-diagonal parts
  n <- nrow(delta.mat)
  K <- length(unique(a))
  total <- rep(0,K+1)
  size <- rep(0,K+1)
  deltamean <- rep(0,K+1)
  names(deltamean) <- c(0,sort(unique(a)))
  for(k in 1:K){
    a_k <- sort(unique(a))[k]
    total[1+k] <- sum(delta.mat[a==a_k,a==a_k])
    size[1+k] <- sum(a==a_k)^2
    deltamean[1+k] <- total[1+k]/size[1+k]
  }
  deltamean[1] <- (sum(delta.mat) - sum(total[-1]))/(n*n - sum(size[-1]))
  return(deltamean)
}

## Trial = split or relocate

Split <- function(a) {
  n <- length(a)
  ua <- unique(a)
  if(length(ua)==n) {
    return(a)
  } else{
    #ua.pick <- sample(x=ua,size=1)
    #a[names(sample(x=which(a==ua.pick),size=1))] <- max(ua)+1
    a.pick <- sample(x=a,size=1)
    pick.ind <- sample(x=which(a==a.pick),size=1)
    a[pick.ind] <- max(a)+1
    return(a)
  }
}

Relocate <- function(a){
  n <- length(a)
  ua <- unique(a)
  if(length(ua)==1) {
    return(a)
  } else{
    pick.ind <- sample(x=1:n,size=1)
    #a.pick <- a[pick.ind]
    #a[pick.ind] <- sample(x=a[-which(a==a.pick)],size=1)
    a[pick.ind] <- sample(x=a[-pick.ind],size=1)
    return(a)
  }
}


##########################
###### Scatterness #######
##########################


scatter = function(dat,cluster.assign,sil_cut=0.1){
  sd_check = apply(dat, 1, sd)

  if(ncol(dat) == 1|!all(sd_check != 0) ){
    dist.dat = as.matrix(dist(dat))
  }else{
    dist.dat = 1 - cor(t(dat))
  }
  sil = silhouette(cluster.assign, dist=dist.dat,diss=T)

  new.dist = dist.dat
  cluster.assign2 = cluster.assign

  while(min(sil[,3]) < sil_cut){
    cluster.assign2<-cluster.assign2[-which(sil[,3] == min(sil[,3]))]
    for (i in unique(cluster.assign2)){
      if(length(which(cluster.assign2==i))==1){
        cluster.assign2<-cluster.assign2[-which(cluster.assign2==i)]
      }
    }
    temp<-cluster.assign2
    for(d in 1:length(cluster.assign2)){
      cluster.assign2[d]<-rank(unique(temp))[which(unique(temp)==temp[d])]
    }#rename cluster index, so it is integer from 1 to k

    new.dist<-new.dist[rownames(new.dist)%in%names(cluster.assign2),
                       colnames(new.dist)%in%names(cluster.assign2)]
    sil <- silhouette(cluster.assign2, dist=new.dist, diss=T)#recalculate silhoutte
  }

  scatter.index = which(!names(cluster.assign)%in%names(cluster.assign2))
  return(scatter.index)
}

##########################
###### Text mining #######
##########################
TextMine <- function(hashtb, pathways, pathway, result, scatter.index=NULL,permutation="all"){
  cat("Performing Text Mining Analysis...\n")
  hashtbAll = hashtb
  k <- length(unique(result))
  if(!is.null(scatter.index)){
    cluster = result
    cluster[scatter.index] = k+1
    hashtb = hashtb[hashtb [,3]%in%which(pathways %in% pathway),]
    tmk = list()
    nperm = 1000
    if (nrow(hashtb) == 0){
      for (i in 1:(k-1)){
        tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)
      }
    }
    else{
      for (i in 1:k){
        e = cluster[cluster == i]
        e = which(pathways %in% names(e))
        hashcl = hashtb[hashtb [,3]%in%e,]
        hashcl = hashcl[duplicated(hashcl[,2]) | duplicated(hashcl[,2], fromLast=TRUE),]
        if (nrow(hashcl) != 0){
          hashf = hashcl
          hashf[,2] = 1
          hashf = aggregate(hashf[,c("row","value")],by = hashf["phrase"],FUN=sum)

          rownames(hashf) = hashf[,"phrase"]
          hashf = hashf[,-1]
          colnames(hashf) = c("count","sum")
          hashap = hashcl
          hashap[,c(2,3,4)] = 0
          mperm = matrix(nrow = nrow(hashf),ncol = nperm)
          for (j in 1:nperm){
            if (permutation=="all"){
              subtb = hashtbAll[hashtbAll[,3]%in%sample(1:length(pathways),length(e)),]
            }else if(permutation=="enriched"){
              subtb = hashtb[hashtb[,3]%in%sample(unique(hashtb[,3]),length(e)),]
            }else{
              stop("Permutation should be 'all' or 'enriched' ")
            }
            subtb = rbind(subtb,hashap)
            subtb = subtb[subtb$phrase %in% hashap$phrase,]
            subtb[,2] = 1
            subtb = aggregate(subtb[,c("row","value")],by = subtb["phrase"],FUN=sum)
            rownames(subtb) = subtb[,"phrase"]
            subtb = subtb[,-1]
            colnames(subtb) = c("count","sum")
            mperm[,j] = subtb[,2]
          }
          hashf[,"p-value"] = apply(cbind(hashf[,2],mperm),1,
                                    function(x)((nperm + 2)-rank(x)[1])/(nperm + 1))
          hashf[,"q-vlaue"] = p.adjust(hashf[,"p-value"],method = "BH")
          tmk[[i]] = hashf[order(hashf[,3],-hashf[,2]),]
        }
        else {tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)}
      }
      tmk[[k+1]] = matrix(NA,nrow = 1,ncol = 4)
    }
  }else{
    hashtb = hashtb[hashtb [,2]%in%which(pathways %in% pathway),]
    tmk = list()
    nperm = 1000
    if (nrow(hashtb) == 0){
      for (i in 1:(k-1)){
        tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)
      }
    }
    else{
      for (i in 1:k){
        e = cluster[cluster == i]
        e = which(pathways %in% names(e))
        hashcl = hashtb[hashtb [,3]%in%e,]
        hashcl = hashcl[duplicated(hashcl[,2]) | duplicated(hashcl[,2], fromLast=TRUE),]
        if (nrow(hashcl) != 0){
          hashf = hashcl
          hashf[,2] = 1
          hashf = aggregate(hashf[,c("row","value")],by = hashf["phrase"],FUN=sum)

          rownames(hashf) = hashf[,"phrase"]
          hashf = hashf[,-1]
          colnames(hashf) = c("count","sum")
          hashap = hashcl
          hashap[,c(2,3,4)] = 0
          mperm = matrix(nrow = nrow(hashf),ncol = nperm)
          for (j in 1:nperm){
            if (permutation=="all"){
              subtb = hashtbAll[hashtbAll[,3]%in%sample(1:length(pathways),length(e)),]
            }else if(permutation=="enriched"){
              subtb = hashtb[hashtb[,3]%in%sample(unique(hashtb[,3]),length(e)),]
            }else{
              stop("Permutation should be 'all' or 'enriched' ")
            }
            subtb = rbind(subtb,hashap)
            subtb = subtb[subtb$phrase %in% hashap$phrase,]
            subtb[,2] = 1
            subtb = aggregate(subtb[,c("row","value")],by = subtb["phrase"],FUN=sum)
            rownames(subtb) = subtb[,"phrase"]
            subtb = subtb[,-1]
            colnames(subtb) = c("count","sum")

            mperm[,j] = subtb[,2]
          }
          hashf[,"p-value"] = apply(cbind(hashf[,2],mperm),1,
                                    function(x)((nperm + 2)-rank(x)[1])/(nperm + 1))
          hashf[,"q-vlaue"] = p.adjust(hashf[,"p-value"],method = "BH")
          tmk[[i]] = hashf[order(hashf[,3],-hashf[,2]),]
        }
        else {tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)}
      }
    }
  }
  return(tmk)
} # End of Text Mining


writeTextOut <- function(tm_filtered,k,pathway.summary) {
  cat("Cluster 1\n", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  write.table(t(rownames(tm_filtered[[1]])[1:15]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
              append = T, row.names=F,col.names=F,na="")
  cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  write.table(t(tm_filtered[[1]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
              append = T, row.names=F,col.names=F,na="")
  cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  write.table(t(tm_filtered[[1]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
              append = T, row.names=F,col.names=F,na="")
  write.table(pathway.summary[[1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T,
              append = T, row.names=F,col.names=F)
for (i in 2:k){
    cat(paste("\nCluster ", i, "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
    cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(rownames(tm_filtered[[i]])[1:15]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
                append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[i]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
                append = T, row.names=F,col.names=F,na="")
    cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[i]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
                append = T, row.names=F,col.names=F,na="")
    write.table(pathway.summary[[i]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T,
                append = T, row.names=F,col.names=F)

  }
}

writeTextOut <- function(tm_filtered,k,pathway.summary,scatter.index=NULL) {
  if(is.null(dim(tm_filtered[[1]]))==TRUE|dim(tm_filtered[[1]])[1] == 0){
    print(paste("No phrase pass q-value threshold in cluster 1"))
    cat("Cluster 1\n", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(pathway.summary[[1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
  } else {
    cat("Cluster 1\n", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(as.character(rownames(tm_filtered[[1]])[1:15])), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[1]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[1]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    write.table(pathway.summary[[1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)

  }

  for (i in 2:k){
    if(is.null(dim(tm_filtered[[i]]))==TRUE|dim(tm_filtered[[i]])[1] == 0){
      print(paste("No phrase pass q-value threshold in cluster ",k,sep = ""))
      cat(paste("\nCluster ", i, "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
      write.table(pathway.summary[[i]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
    } else {
      cat(paste("\nCluster ", i, "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
      cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
      write.table(t(as.character(rownames(tm_filtered[[i]])[1:15])), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
      write.table(t(tm_filtered[[i]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
      write.table(t(tm_filtered[[i]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(pathway.summary[[i]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
    }

  }
  if(!is.null(scatter.index)){
    cat(paste("\nSingleton Term", "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
    write.table(pathway.summary[[k+1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T,
                append = T, row.names=F,col.names=F)
  }
}

textMine <- function(hashtb,pathways,cluster.assign,scatter.index=NULL,thres=0.05,permutation="all"){
  tmk <- TextMine(hashtb=hashtb, pathways= pathways,
                  pathway=names(cluster.assign), result=cluster.assign, scatter.index,permutation=permutation)
  C <- length(unique(cluster.assign))
  tm_filtered <- list()
  for (i in 1:C){
    tm_filtered[[i]] <- tmk[[i]][which((as.numeric(tmk[[i]][,4]) < thres)), ]
  }
  if(!is.null(scatter.index)){
    tm_filtered[[C+1]] = tmk[[C+1]]
    cluster = cluster.assign
    cluster[scatter.index] = C+1
    pathway.summary <- lapply(1:(C+1), function(x) names(which(cluster==x)))
    writeTextOut(tm_filtered,C,pathway.summary,scatter.index=scatter.index)
  }else{
    pathway.summary <- lapply(1:C, function(x) names(which(cluster.assign==x)))
    writeTextOut(tm_filtered,C,pathway.summary)
  }
  return(tm_filtered)
}

##########################
###### ACS/ADS_DE plot ###
##########################
ARS_to_size <- function(ARSp,factor=2){
  if(ARSp > 0.05) {
    return(1)
  } else {
    return(-log10(ARSp)*factor)
  }
}

ACS_ADS_DE <- function(ds1,ds2,DEevid1,DEevid2,ACSp,ADSp,cluster=NULL,
                       highlight.pathways=NULL,lb=0,ub=1,size.scale=4){

  P <- length(ACSp)
  ACS_size=sapply(ACSp,ARS_to_size)
  ADS_size=sapply(ADSp,ARS_to_size)


  if(!is.null(cluster)){
    RB = rainbow(length(unique(cluster)))
    color = c()
    for (i in 1:length(cluster)) {
      if(cluster[i] != "scatter"){
        color[i] = RB[as.numeric(cluster[i])]
      }else{
        color[i] = "grey50"
      }
    }
  }else if(!is.null(highlight.pathways)){
    color = rep("grey50",P)
    color[highlight.pathways] = "red"
  }else{
    color = rep("black",P)
  }


  if(!is.null(highlight.pathways)){
    index = 1:P
    index[-highlight.pathways] = ""
  }else{
    index = rep("",P)
  }
  data_ACS <- data.frame(ds1_score=DEevid1,ds2_score=DEevid2,
                         ACS_size=ACS_size,index=index,
                         color_pos=color)
  data_ADS <- data.frame(ds1_score=DEevid1,ds2_score=DEevid2,
                         ADS_size=ADS_size,index=index,
                         color_neg=color)

  p_pos <-ggplot(data_ACS, aes(x=ds2_score, y=ds1_score,label=index)) +
    geom_point(size = data_ACS$ACS_size*size.scale, shape=16, color = data_ACS$color_pos)+
    geom_text(size=8*size.scale,parse=TRUE,color="black",hjust = -0.05,vjust=-0.05) +
    theme_bw() +
    coord_fixed(ylim=c(lb,ub),xlim=c(lb,ub)) +
    labs(x="",y="") +
    scale_x_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) +
    scale_y_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) +
    theme(#legend.title = element_blank(),
      axis.line = element_line(colour = "black"),
      #axis.line=element_blank(),
      axis.text.x = element_text(size = 40,face = "bold"),
      axis.text.y = element_text(size = 40,face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_line(linetype = 'solid',#size = 2,
                                      colour = "white"),
      panel.grid.minor = element_line(linetype = 'solid',
                                      colour = "white"),
      panel.background = element_rect(fill = "#FFF1E1")) +
    annotate("text", x = (ub-0.15), y = lb, fontface=2,
             label = paste(ds2,sep=""),
             size=8*size.scale,colour="blue",hjust=0.6,vjust=0.1) +
    annotate("text", x = lb, y = (ub-0.15), fontface=2,
             label=paste(ds1,sep=""),
             size=8*size.scale,colour="blue",vjust=0,hjust=0.2)

  p_neg <-ggplot(data_ADS, aes(x=ds1_score, y=ds2_score,label=index)) + #label=index
    geom_point(size = data_ADS$ADS_size*size.scale, shape=16, color = data_ADS$color_neg)+
    geom_text(size=8*size.scale,parse=TRUE,color="black",hjust = -0.05,vjust=-0.05) +
    theme_bw() +
    coord_fixed(ylim=c(lb,ub),xlim=c(lb,ub)) +
    labs(x="",y="") +
    scale_x_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) +
    scale_y_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) +
    theme(#legend.title = element_blank(),
      axis.line = element_line(colour = "black"),
      #axis.line=element_blank(),
      axis.text.x = element_text(size = 40,face = "bold"),
      axis.text.y = element_text(size = 40,face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_line(linetype = 'solid',#size = 2,
                                      colour = "white"),
      panel.grid.minor = element_line(linetype = 'solid',
                                      colour = "white"),
      panel.background = element_rect(fill = "#EBF5FF")) +
    annotate("text", x = (ub-0.15), y = lb, fontface=2,
             label = paste(ds1,sep=""),
             size=8*size.scale,colour="blue",hjust=0.6,vjust=0.1) +
    annotate("text", x = lb, y = (ub-0.15), fontface=2,
             label=paste(ds2,sep=""),
             size=8*size.scale,colour="blue",vjust=0,hjust=0.2)

  #ggsave(filename=paste(ds2,"_",ds1,"_ACS_figure",".pdf",sep=""),p_pos,
  #width = 10, height = 10)

  #ggsave(filename=paste(ds1,"_",ds2,"_ADS_figure",".pdf",sep=""),p_neg,
  #width = 10, height = 10)

  plist <- list(p_pos,p_neg)
  return(plist)
}

##########################
######    parseXML     ###
##########################
parseRelation <- function(pathwayID, keggSpecies="hsa", binary = T, sep = "-") {
  # download xml file
  download.kegg(pathway.id = pathwayID, keggSpecies, kegg.dir = ".", file.type="xml")
  # generate relation matrix
  pathName = unlist(as.list(KEGGPATHID2NAME)[pathwayID])

  xmlFile = paste0(getwd(), "/",keggSpecies, pathwayID,".xml")
  pathway = KEGGgraph::parseKGML(xmlFile)
  pathway = KEGGgraph::splitKEGGgroup(pathway)

  entries = KEGGgraph::nodes(pathway)
  types = sapply(entries, getType)
  relations = unique(KEGGgraph::edges(pathway)) ## to avoid duplicated edges
  relationNum = length(relations)
  entryNames = as.list(sapply(entries, getName))
  if(any(types == "group") || any(types=="map")){
    entryNames = entryNames[!(types %in% c("group","map"))]
  }
  entryIds = names(entryNames)
  entryNames = lapply(1:length(entryNames), function(i) paste(entryNames[[i]],collapse=sep))
  names(entryNames) = entryIds

  entryNames.unique = unique(entryNames)
  entryNum = length(entryNames.unique)

  relation.mat = matrix(0, entryNum, entryNum)
  rownames(relation.mat) = colnames(relation.mat) = entryNames.unique

  ## if no relation edge, just return
  if(relationNum == 0){
    print(paste0("There is no topological connected gene nodes in ", pathName))
    return(relation.mat)
  }

  entry1 = getEntryID(relations)[,1]
  entry2 = getEntryID(relations)[,2]
  for(i in 1:length(relations)){
    if(entry1[i] %in% names(entryNames) && entry2[i] %in% names(entryNames)){
      relation.mat[entryNames[[entry1[i]]],entryNames[[entry2[i]]]]=1
    }
    else{
      print(paste("connections not included:",entry1[i], entry2[i], sep=" "))
    }
  }

  file.remove(xmlFile)
  return(relation.mat)
}
##########################
###### KEGG module SA ####
##########################
SA_module_M = function(sp.mat, xmlG, M, nodes, B = 1000,
                       G.ini.list=NULL, reps_eachM = 100,topG_from_previous=10,
                       Tm0=10,mu=0.95,epsilon=1e-5,
                       N=1000,run=10000,seed=12345,sub.num=1){
  #Null distribution for M
  set.seed(seed)
  null.sp.dist = rep(NA,B)
  for(b in 1:B){
    permute.set <- sample(xmlG,M)
    permute.mat <- sp.mat[match(permute.set,row.names(sp.mat)),
                          match(permute.set,colnames(sp.mat))]
    null.sp.dist[b] <- mean(c(permute.mat[lower.tri(permute.mat)]))
  }
  null.sp.mean = mean(null.sp.dist)
  null.sp.median = median(null.sp.dist)
  if(is.null(G.ini.list)){
    p.mean.ls = c()
    G.module.ls = list()
    SP.ls = c()
    for (l in 1:reps_eachM) {
      ##Initialize
      if(length(nodes) == M){
        G.module = nodes
      }else{
        G.module = sample(nodes, M)
      }
      SPc = avgSP(G.module, sp.mat)
      r = 0
      count = 0
      Tm = Tm0
      while((length(nodes)>M) & (r < run) & (count < N) & (Tm >= epsilon)) {
        #pi = exp(-GPc/Tm) ## Boltzmann dist #may need a different Tm or -logP to be comparable?
        #print(SPc)
        ##New trial
        r = r+1
        a.node = sample(setdiff(nodes,G.module),sub.num)
        G.module_new = G.module
        G.module_new[sample(M,sub.num)] = a.node

        SPn = avgSP(G.module_new, sp.mat)

        if(SPn < SPc | SPc == Inf) {
          ##accept
          SPc = SPn;
          G.module = G.module_new;
        }else{
          count = count + 1;
          p = exp((SPc-SPn)/Tm)
          #print(p)
          r = min(1,p); ## acceptance prob.
          u <- runif(1);
          if(u>r) {
            ##not accept
            Tm <- Tm*mu
          } else {
            SPc = SPn;
            G.module = G.module_new;
          }
        }
      }
      G.module.ls[[l]] = G.module
      p.mean.ls[l] = (sum(null.sp.dist<= SPc) + 1)/(B+1)
      SP.ls[l] = SPc
    }
  }else{
    each.times = round(reps_eachM/topG_from_previous)
    case.index = expand.grid(1:length(G.ini.list),1:each.times)
    p.mean.ls = c()
    G.module.ls = list()
    SP.ls = c()
    for (l in 1:nrow(case.index)) {
      G.ini = G.ini.list[[case.index[l,1]]]
      ##Initialize
      if(length(nodes) == M){
        G.module = nodes
      }else{
        x = M-length(G.ini)
        G.module = c(G.ini,sample(setdiff(nodes,G.ini),x))
      }
      SPc = avgSP(G.module, sp.mat)
      r = 0
      count = 0
      Tm = Tm0
      while((length(nodes)>M) & (r < run) & (count < N) & (Tm >= epsilon)) {
        #pi = exp(-GPc/Tm) ## Boltzmann dist #may need a different Tm or -logP to be comparable?
        #print(SPc)
        ##New trial
        r = r+1
        a.node = sample(setdiff(nodes,G.module),sub.num)
        G.module_new = G.module
        G.module_new[sample(M,sub.num)] = a.node

        SPn = avgSP(G.module_new, sp.mat)

        if(SPn < SPc | SPc == Inf) {
          ##accept
          SPc = SPn;
          G.module = G.module_new;
        }else{
          count = count + 1;
          p = exp((SPc-SPn)/Tm)
          #print(p)
          r = min(1,p); ## acceptance prob.
          u <- runif(1);
          if(u>r) {
            ##not accept
            Tm <- Tm*mu
          } else {
            SPc = SPn;
            G.module = G.module_new;
          }
        }
      }
      G.module.ls[[l]] = G.module
      p.mean.ls[l] = (sum(null.sp.dist<= SPc) + 1)/(B+1)
      SP.ls[l] = SPc
    }

  }
  p.sd.ls = sqrt(p.mean.ls*(1-p.mean.ls)/B)
  r.p = rank(p.mean.ls,ties.method = "random")
  index = match(1:topG_from_previous,r.p)
  best.index = which(r.p == 1)

  minG = G.module.ls[[best.index]]
  sp = SP.ls[best.index]
  p.mean = p.mean.ls[best.index]
  p.sd = p.sd.ls[best.index]

  top.G = G.module.ls[index]
  top.pmean = p.mean.ls[index]
  top.psd = p.sd.ls[index]
  top.sp = SP.ls[index]

  return(list(minG = minG,sp = sp,p.mean = p.mean,p.sd = p.sd,
              top.G = top.G,top.sp = top.sp,top.pmean = top.pmean,top.psd = top.psd,
              null.sp.mean = null.sp.mean,null.sp.median = null.sp.median))
}

avgSP = function(G.module, sp.mat){
  m = length(G.module)
  set.mat = sp.mat[match(G.module,row.names(sp.mat)),
                   match(G.module,colnames(sp.mat))]
  G.sp = mean(c(set.mat[lower.tri(set.mat)]))
  return(G.sp)
}
