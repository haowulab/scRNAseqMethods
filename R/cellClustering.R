##############################################
## Functions for scRNA-seq clustering
## Inputs are **normalized counts** in the raw scale.
## Features are not selected
##
## ouput are clusters
##############################################


#' Wrapper function for scRNA-seq cell clustering
#' @param X A matrix of *normalized* counts (in raw scale), rows are genes, columns are cells
#' @param method Cell clustering method to be used
#' @param K number of clusters, an integer
#' @return A vector, with length being the number of cells, for the cluster assignment. 
#' @note The input X should contain all genes. Features are not selected.
#' @export

cellClustering <- function(x,
                           method = c("SC3", "Seurat", "monocle", "cidr", 
                                      "TSCAN", "SHARP"),
                           K ) {
    method = match.arg(method)
    switch(method,
           SC3 = cluster.SC3(x, K),
           Seurat = cluster.Seurat(x, K),
           monocle = cluster.monocle(x, K),
           cidr = cluster.cidr(x, K),
           TSCAN = cluster.TSCAN(x, K),
           SHARP = cluster.SHARP(x, K)
           )
}

### clustering using Seurat
cluster.Seurat <- function(counts, K) {
    require(Seurat)
    ## PCA dimension reduction
    seuset = CreateSeuratObject( counts )
    seuset = NormalizeData(object = seuset)
    seuset = FindVariableFeatures(object = seuset)
    seuset = ScaleData(object = seuset)
    seuset = RunPCA(object = seuset)
    seuset = FindNeighbors(object = seuset)
    seuset = FindClusters(object = seuset)
    return(seuset@active.ident)
}

### clustering using SC3
cluster.SC3 <- function(counts, K) {
    require(SC3)
    sce = SingleCellExperiment(
        assays = list(
            counts = as.matrix(counts),
            logcounts = log2(as.matrix(counts) + 1)
        )
    )
    rowData(sce)$feature_symbol <- rownames(sce)
    sce = sce[!duplicated(rowData(sce)$feature_symbol), ]
    sce = sc3_prepare(sce)
    if( missing(K) ) { ## estimate number of clusters
        sce = sc3_estimate_k(sce)
        K = metadata(sce)$sc3$k_estimation
    }

    sce = sc3_calc_dists(sce)
    sce = sc3_calc_transfs(sce)
    sce = sc3_kmeans(sce, ks = K)
    sce = sc3_calc_consens(sce)
    colTb = colData(sce)[,1]
    return(colTb)
}


### Monocle
cluster.monocle <- function(counts, K) {
    require(monocle)
    dat <- newCellDataSet(counts)
    dat <- estimateSizeFactors(dat)
    dat <- estimateDispersions(dat)

    ## pick marker genes for cell clustering
    disp_table <- dispersionTable(dat)
    ix <- order(disp_table[,"mean_expression"], decreasing=TRUE)
    unsup_clustering_genes <- disp_table[ix[50:1000], "gene_id"]
    ## unsup_clustering_genes <- subset(disp_table, mean_expression >= 1)
    dat <- setOrderingFilter(dat, unsup_clustering_genes)
    ## the follwoing step can be slow. Need to keep marker genes number low.
    dat <- reduceDimension(dat, reduction_method="tSNE")

    ## clustering
    if( !missing(K) )
        dat <- clusterCells(dat, num_clusters = K)
    else
        dat <- clusterCells(dat)

    pData(dat)$Cluster
}

### CIDR
cluster.cidr <- function(counts, K) {
    require(cidr)
    sData <- scDataConstructor(counts)
    sData <- determineDropoutCandidates(sData)
    sData <- wThreshold(sData)
    sData <- scDissim(sData)
    sData <- scPCA(sData, plotPC=FALSE)
    sData <- nPC(sData)
    if(missing(K))
        sData <- scCluster(sData)
    else
        sData <- scCluster(sData, nCluster=K)

    return(sData@clusters)
}

### TSCAN
cluster.TSCAN <- function(counts, K) {
    require(TSCAN)
    procdata <- preprocess(counts, cvcutoff=0.5)
    if( !missing(K) )
        lpsmclust <- exprmclust(procdata, clusternum=K)
    else
        lpsmclust <- exprmclust(procdata)

    return(lpsmclust$clusterid)
}

### SHARP
cluster.SHARP <- function(counts, K) {
    require(SHARP)
    data = log2(counts+1)
    clusters.SHARP = SHARP(data, N.cluster=K)
    return(clusters.SHARP[[1]])
}

