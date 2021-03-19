#############################################
## Functions for DE
#############################################

#' Wrapper function for scRNA-seq DE
#' @param X A matrix of *raw* counts, rows are genes, columns are cells
#' @param method DE method to be used
#' @param ... Other parameters for the corresponding method.
#' @return An imputed matrix. This could be in raw or log scales, 
#' depending on the method. 
#' @export

DE <- function(X, method=c("MAST", "SC2P"), ...) {
  method = match.arg(method)
  switch(method,
         DrImpute = DE.DrImpute(X, ...),
         SAVER = DE.SAVER(X, ...)
  )
}

DE.DrImpute <- function(x, ...) {
  
}

DE.SAVER <- function(X, ...) {
  
}


