#############################################
## Functions for simulator
#############################################

#' Wrapper function for scRNA-seq simulator
#' @description Given a template dataset, simulate scRNA-seq counts. 
#' The template data must be from one specific cell types. 
#' @param X A matrix of *raw* counts for one cell type. Rows are genes, columns are cells.
#' @param ncells Number of cells in the simulated data
#' @param totalCounts Total read counts for cells. Can be an integer 
#' (then all cells will have the same counts), or a vector of integers with 
#' length being the number of cells.
#' @param method simulator method to be used
#' @param ... Other parameters for the corresponding method.
#' @return A simulated count matrix. Rows are genes, columns are cells. 
#' Number of cells is specified by the user. 
#' Number of genes is the same as the input.
#' @export

simulator <- function(X, ncells, totalCounts, 
                      method=c("scDesign", "splat", "POWSC"), ...) {
  method = match.arg(method)
  switch(method,
         scDesign = simulator.scDesign(X, ncells, totalCounts, ...),
         splat = simulator.splat(X, ncells, totalCounts, ...),
         POWSC = simulator.POWSC(X, ncells, totalCounts, ...)
  )
}

simulator.POWSC <- function(X, ncell, ...) {
  require(POWSC)
  est_paras = Est2Phase(sce = X)
  POWSC_sim = Simulate2SCE(n = ncell, estParas1 = est_paras, estParas2 = est_paras)
  return(assays(POWSC_sim$sce)$counts)
}

simulator.scDesign <- function(X, ncells, ...) {
  require(scDesign)
  scDesign_sim = design_data(realcount = X, S = 1e7, 
                             ncell = ncells, ngroup = 1) 
  return(scDesign_sim) 
}

simulator.splat <- function(X, ncells, ...) {
  require(splatter)
  splat_paras = splatEstimate(X)
  splat_paras@nCells = ncells
  splat_paras@batchCells = ncells
  splat_sim = splatSimulate(splat_paras, method = "single")
  return(assays(splat_sim)$counts)
}


