# scRNAseqMethods
An R package with some wrapper functions for analyzing scRNA-seq data. 

The package contains functions for following tasks
- normalization
- batch effect removal
- imputation
- unsupervised cell clustering
- supervised cell typing
- data simulator

# Installation
## Installing scRNAseqMethods package 
```{r, message=FALSE, eval=FALSE}
library(devtools)
install_github("haowulab/scRNAseqMethods", build_vignettes=TRUE)
```

## Installing other packages 
`scRNAseqMethods` provides wrapper functions for the methods implemented in 
many other packages. The scRNAseqMethods package itself doesn't depend on 
those packages, but you need to have specific packages installed if you
want to use their methods. 

We provide a script to install all those packages. It could take a little while
since there might be other dependence packages need to be installed. 
Also there could be installation errors. In that case, please take a look
at the `installPkg` function and run it line by line manually. 

```{r, eval=FALSE}
library(scRNAseqMethods)
installPkg()
```

# Examples
There is a small dataset distributed with the package for illustrative purpose.
It is a human pancreas data set with 3 subjects, each has 6 cell types. 
We did pre-filtering and kept a total of 2000 genes and 702 cells. 

```{r, message=FALSE, eval=FALSE}
library(scRNAseqMethods)
data(Seg_PancreasData)
```

## Data normalization 
The `normalization` function takes a raw count matrix 
(rows are genes, columns are cells) and a method, 
and returns a normalized count matrix of the same dimension. 

```{r, eval=FALSE}
Xnorm = normalization(Seg_counts, method="total")
Xnorm = normalization(Seg_counts, method="scran")
```

## Batch effect correction
The `batchCorrection` function takes a raw count matrix 
(rows are genes, columns are cells), a factor vector, and a method. 
It returns a matrix of the same dimension with batch effect corrected. 
Note that the return matrix could be in raw or log scales, 
depending on the method. 

```{r, eval=FALSE}
Xcorrect = batchCorrection(Seg_counts, Seg_subject_ID, method="harmony")
```
## Imputation
The `imputaton` function takes a raw count matrix, a method, 
and return a matrix of the same dimension with missing data imputed. 
```{r, eval=FALSE}
Ximpute = imputation(Seg_counts, method="DrImpute")
```

## Cell clustering 
The `cellClustering` function takes a raw count matrix, a method, and 
number of desired clusters, then return a vector for clustering results.
```{r, eval=FALSE}
res = cellClustering(Seg_counts, method="SHARP", K=6)
table(res, Seg_true_cell_label)
```

## Cell typing


## Data simulator
The `simulator` function take a template data (raw counts for **one** cell type),
number of cells, total counts for each cell, and generate a 
simulated count matrix. 

```{r, eval=FALSE}
## take the counts for alpha cells in subject 1
ix = Seg_subject_ID == 1
Y0 = Seg_counts[,ix]
ix.celltype = Seg_true_cell_label[ix]=="alpha"
mat = Y0[,ix.celltype]
ncell = 100
Xsim = simulator(mat, ncell, totalCounts=5000, method="POWSC")
dim(Xsim)
## check the correlation between simulated and template data
cor(mat[,1:5], Xsim[,1:5])
```

