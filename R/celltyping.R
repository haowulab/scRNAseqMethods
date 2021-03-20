#############################################
## Functions for celltyping
#############################################

#' Wrapper function for scRNA-seq celltyping
#' @param X A matrix of *raw* counts rows are genes, columns are cells
#' @param method celltyping method to be used
#' @param ... Other parameters for the corresponding method.
#' @return An imputed matrix. This could be in raw or log scales, 
#' depending on the method. 
#' @export

celltyping <- function(Xtrain, labelsTraing, Xtest, 
                       method=c("SingleR", "scmapCluster", "scmapCell",
                                "CHETAH", "Garnett"), ...) {
  method = match.arg(method)
  switch(method,
         singleR = celltyping.SingleR(Xtrain, labelsTraing, Xtest, ...),
         singleR = celltyping.scmapCluster(Xtrain, labelsTraing, Xtest, ...),
  )
}

## singleR
celltyping.SingleR <- function(Xtrain, labelsTraing, Xtest, ...) {
  require(SingleR)
  pred.hesc <- SingleR(test = Xtest, ref = Xtrain, 
                       labels = labelsTraing)
  return(pred.hesc$labels)
}

celltyping.scmapCluster <- function(Xtrain, labelsTraing, Xtest, ...) {
  require(scmap)
  
}

celltyping.CHETAH <- function(Xtrain, labelsTraing, Xtest, ...) {
  require(CHETAH)
}


