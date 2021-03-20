###########################################################
## A function to install all packages used in this pacakge 
###########################################################

#' Install all packages used in this package
#' @description The packages used include 
#' For data normalization: DESeq2, scran, SC2P, SCnorm 
#' For batch effect correction: batchelor (MNN), harmony
#' For imputation: SAVER, DrImpute
#' For cell clustering: SC3, Seurat, monocle3, cidr, TSCAN, SHARP
#' For differential expression: MAST, SC2P
#' For cell typing: 
#' For simulator: POWSC, splatter, scDesign
#' @export

installPkg <- function() {
  require(devtools)
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)

  BiocManager::install("SingleCellExperiment")
  BiocManager::install("scater")
  
  BiocManager::install("SCnorm")
  BiocManager::install("scran")
  BiocManager::install("DESeq2")
  
  BiocManager::install("batchelor")
  install_github("immunogenomics/harmony", 
                 build_vignettes = T, dependencies = T)
  
  BiocManager::install("Seurat")
  BiocManager::install("SC3")
  devtools::install_github("VCCRI/CIDR", build_vignettes = T, dependencies = T)
  BiocManager::install("TSCAN")
  install_github("shibiaowan/SHARP", build_vignettes = T, dependencies = T)
  devtools::install_github('cole-trapnell-lab/monocle3',
                           build_vignettes = T, dependencies = T)
  
  install.packages("SAVER")
  install.packages("DrImpute")
  
  BiocManager::install("MAST")
  BiocManager::install("splatter")
  install_github("Vivianstats/scDesign", build_vignettes = T, dependencies = T)
  install_github("haowulab/SC2P", build_vignettes=TRUE) 
  # POWSC requires SC2P to be installed first
  install_github("suke18/POWSC", build_vignettes = T, dependencies = T)
  
}
