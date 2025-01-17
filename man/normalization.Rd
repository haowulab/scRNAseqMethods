% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{normalization}
\alias{normalization}
\title{Wrapper function for scRNA-seq data normalization}
\usage{
normalization(X, method = c("total", "DESeq", "scran", "SC2P", "SCnorm"))
}
\arguments{
\item{X}{A matrix of raw counts, rows are genes, columns are cells}

\item{method}{Normalization method to be used}
}
\value{
A normalized data matrix, in raw scale.
}
\description{
Wrapper function for scRNA-seq data normalization
}
