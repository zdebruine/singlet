#' Load the pbmc3k dataset
#'
#' 2,700 peripheral blood mononuclear cells (PBMC) from 10x genomics taken from the "SeuratData" package
#'
#' @description
#' This dataset is adapted directly from the Satija lab "pbmc3k" dataset used in their popular tutorial on guided clustering. It is provided in this package for convenience since "SeuratData" is not available on CRAN.
#'
#' For more information, please see their documentation.
#'
#' @returns Seurat object with \code{$cell_type} info in the \code{meta.data} slot.
#'
#' @export
#'
get_pbmc3k_data <- function() {
  data(pbmc3k)
  pbmc3k
  A <- CreateSeuratObject(counts = new("dgCMatrix", i = pbmc3k$i, p = pbmc3k$p, Dim = pbmc3k$Dim, Dimnames = pbmc3k$Dimnames, x = as.numeric(inverse.rle(pbmc3k$x))))
  A@meta.data$cell_type <- pbmc3k$cell_type
  A
}
