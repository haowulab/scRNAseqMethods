###########################################################
## scRNA-seq normalization methods
###########################################################

#' Wrapper function for scRNA-seq data normalization
#' @param X A matrix of raw counts, rows are genes, columns are cells
#' @param method Normalization method to be used
#' @return A normalized data matrix, in raw scale. 
#' @export

normalization <- function(X, 
                          method=c("total", "DESeq", "scran", "SC2P", "SCnorm")) {
    method = match.arg(method)
    res <- switch(method,
                  total = normalizeSC.total(X),
                  DESeq = normalizeSC.DESeq(X),
                  scran = normalizeSC.scran(X),
                  SC2P = normalizeSC.SC2P(X),
                  SCnorm = normalizeSC.SCnorm(X))
    res
}


##### Use total counts to normalize
normalizeSC.total <- function(X) {
    k=colSums(X)
    k = k/median(k)
    Xnorm = sweep(X, 2, k, FUN="/")
    return(Xnorm)
}

##### Use MR in DESeq
normalizeSC.DESeq <- function(X) {
    require(DESeq2)
    k = estimateSizeFactorsForMatrix(X)
    Xnorm = sweep(X, 2, k, FUN="/")
    return(Xnorm)
}


##### use SC2P
normalizeSC.SC2P <- function(X) {
    require(SC2P)
    ## make a dummy design
    design = rep(1, ncol(X))
    colnames(X) <- names(design) <- paste0("sample", 1:ncol(X))
    phenoData <- new("AnnotatedDataFrame", data=data.frame(design))
    eset <- ExpressionSet(assayData=X, phenoData=phenoData)
    normdata <- eset2Phase(eset)
    offset <- 2^assayData(normdata)$Offset
    Xnorm <- (assayData(normdata)$exprs+1) / offset
    return(Xnorm)
}

##### use SCnorm - this is crazily slow
normalizeSC.SCnorm <- function(X) {
    require(SCnorm)
    require(parallel)
    ncores = detectCores()-3
    if(ncores < 1) ncores = 1

    ## make a dummy design
    design = rep(1, ncol(X))
    DataNorm <- SCnorm(Data = X,
                       Conditions = design,
                       PrintProgressPlots = TRUE,
                       FilterCellNum = 10,
                       NCores=ncores)
    Xnorm <- normcounts(DataNorm)
    return(Xnorm)

}


##### Use scran -- this is the John C. Marioni 2016 GB paper
normalizeSC.scran <- function(X) {
    require(scran)
    k <- calculateSumFactors(X)
    Xnorm = sweep(X, 2, k, FUN="/")
    return(Xnorm)
}

