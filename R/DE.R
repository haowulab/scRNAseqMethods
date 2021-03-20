#############################################
## Functions for DE
#############################################

#' Wrapper function for scRNA-seq DE. This only implements the simple
#' two-group comparison case. We can certainly do others, but the syntax
#' will be complicated. I will keep it simple for now.
#' @param X A matrix of *raw* counts, rows are genes, columns are cells
#' @param design A vector representing the treatment groups.
#' It must be a vector of 0 and 1. The length of the vector must match the
#' number of columns of input count matrix
#' @param ... Other parameters for the corresponding method.
#' @return Depends on the method. Usually a data frame, rows are for genes.
#' Columns are different information including p-values, test statistics, etc.
#' @export


DE <- function(X, design, method=c("MAST", "SC2P"), ...) {
  method = match.arg(method)
  switch(method,
         MAST = DE.MAST(X, design, ...),
         SC2P = DE.SC2P(X, design, ...)
  )
}

DE.MAST <- function(x, design, ...) {
  ## normalize and take log
  Xnorm = normalization(x, method="total")
  ltpm = log(Xnorm + 1)

  sca <- FromMatrix(ltpm, cData=data.frame(design=as.factor(design)))
  cdr2 <- colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  thres <- thresholdSCRNACountMatrix(assay(sca),
                                     nbins=200, min_per_bin=30)
  assays(sca) <- list(thresh=thres$counts_threshold,
                      tpm=assay(sca))

  ## fit model and perform test
  fit <- zlm(~design, sca)
  lrt <- lrTest(fit, "design")

  ## return p-values only
  return(lrt[,,'Pr(>Chisq)'])
}

DE.SC2P <- function(x, design, ...) {
  phenoData <- new("AnnotatedDataFrame",
                   data=data.frame(design=as.factor(design)))
  colnames(x) = 1:length(design)
  eset <- ExpressionSet(assayData=x, phenoData=phenoData)

  ## estimate phases
  data <- eset2Phase(eset)
  de.sc2p <- twoPhaseDE(data, design="design", test.which=1, offset="sf")
  return(de.sc2p)
}


