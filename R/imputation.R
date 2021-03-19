#############################################
## Functions for imputation
#############################################

#' Wrapper function for scRNA-seq imputation
#' @param X A matrix of *raw* counts, rows are genes, columns are cells
#' @param method Imputation method to be used
#' @param ... Other parameters for the corresponding method.
#' @return An imputed matrix. This could be in raw or log scales, 
#' depending on the method. 
#' @export
imputation <- function(X, method=c("DrImpute", "SAVER"), ...) {
  method = match.arg(method)
  switch(method,
         DrImpute = imputation.DrImpute(X, ...),
         SAVER = imputation.SAVER(X, ...)
  )
}

imputation.DrImpute <- function(x, ...) {
  
}

imputation.SAVER <- function(X, ...) {
  
}



