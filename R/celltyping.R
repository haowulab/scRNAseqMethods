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

celltyping <- function(X, method=c("MAST", "SC2P"), ...) {
  method = match.arg(method)
  switch(method,
         DrImpute = celltyping.DrImpute(X, ...),
         SAVER = celltyping.SAVER(X, ...)
  )
}

celltyping.DrImpute <- function(x, ...) {
  
}

celltyping.SAVER <- function(X, ...) {
  
}


