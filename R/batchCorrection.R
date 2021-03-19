#############################################
## Functions for batch effect correction
#############################################

#' Wrapper function for scRNA-seq batch effect correction
#' @param X A matrix of counts (in raw scale), rows are genes, columns are cells
#' @param method Cell clustering method to be used
#' @param batchID A fact, with length being the number of cells, for the batch. 
#' @param ... Other parameters for the corresponding method.
#' @return A batch-corrected matrix. This could be in raw or log scales, depending on the method. 
#' @export
 
batchCorrection <- function(X, method=c("MNN", "harmony"), ...) {
    method = match.arg(method)
    switch(method,
           MNN = batchCorrection.MNN(x, K, ...),
           harmony = batchCorrection.harmony(x, K, ...)
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
 
