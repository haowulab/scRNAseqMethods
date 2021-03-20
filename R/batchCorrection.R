#############################################
## Functions for batch effect correction
#############################################

#' Wrapper function for scRNA-seq batch effect correction
#' @param X A matrix of counts (in raw scale), rows are genes, columns are cells
#' @param batchID A vector of factors, with length being the number of cells, 
#' for the batch. 
#' @param method Cell clustering method to be used
#' @param ... Other parameters for the corresponding method.
#' @return A batch-corrected matrix. This could be in raw or log scales, depending on the method. 
#' @export
 
batchCorrection <- function(X, batchID, method=c("MNN", "harmony"), ...) {
    method = match.arg(method)
    switch(method,
           MNN = batchCorrection.MNN(X, batchID, ...),
           harmony = batchCorrection.harmony(X, batchID, ...)
    )
}


## MNN
batchCorrection.MNN <- function(X, batchID, ...) {
    require(batchelor)
    res = mnnCorrect(X, batchID, ...)
    return(assays(res)$corrected)
}

## harmony
batchCorrection.harmony <- function(X, batchID, ...) {
    require(harmony)
    ## first normalize data to get log normalized counts
    Xnorm = normalization(X, method="total")
    Xcorrected <- HarmonyMatrix(Xnorm, batchID, do_pca = FALSE, ...)
    return(Xcorrected)
}
 
