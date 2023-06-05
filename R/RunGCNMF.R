#' Run Graph-Convolutional Non-negative Matrix Factorization
#'
#' @description Run NMF with weighted convolution determined by edges in a graph of dimensions \code{n x n}, where \code{n} is the number of columns in the matrix.
#'
#' @inheritParams RunNMF
#' @param graph A graph to use, either directed or undirected
#' @param verbose print updates to console
#' @param k rank of the factorization (no automatic rank determination for GCNMF. Use \code{\link{RunNMF}}). Alternatively, specify an initial \code{w} matrix of dimensions \code{m x k}, where \code{m} is the number of rows in the matrix to be factorized.
#'
#' @return Returns a Seurat object with the GCNMF model stored in the reductions slot
#'
#' @details Use \code{set.seed()} to guarantee reproducibility!
#' @rdname RunGCNMF
#' @aliases RunGCNMF.Seurat
#' @name RunGCNMF.Seurat
#'
#' @seealso \code{\link{RunNMF}}
#'
#' @export
#'
RunGCNMF.Seurat <- function(object,
                          graph,
                          k,
                          split.by = NULL,
                          assay = NULL,
                          tol = 1e-5,
                          L1 = 0.01,
                          L2 = 0,
                          verbose = 2,
                          reduction.name = "gcnmf",
                          reduction.key = "GCNMF_",
                          maxit = 100,
                          threads = 0,
                          features = NULL,
                          ...) {
  if (is.null(assay)) {
    assay <- names(object@assays)[[1]]
  }
  
  # check if data has been normalized
  v <- object@assays[[assay]]@data@x
  if (sum(as.integer(v)) == sum(v)) {
    object <- PreprocessData(object, assay = assay)
  }
  A <- object@assays[[assay]]@data
  
  if (!is.null(features)) {
    if (features[[1]] == "var.features") {
      A <- A[object@assays[[assay]]@var.features, ]
    } else if (is.integer(features) || is.character(features)) {
      # array of indices or rownames
      A <- A[features, ]
    } else {
      stop("'features' vector was invalid.")
    }
  }
  
  rnames <- rownames(A)
  cnames <- colnames(A)
  
  if (!is.null(split.by)) {
    split.by <- as.integer(as.numeric(as.factor(object@meta.data[[split.by]]))) - 1
    if (any(sapply(split.by, is.na))) {
      stop("'split.by' cannot contain NA values")
    }
    A <- weight_by_split(A, split.by, length(unique(split.by)))
  }
  At <- Matrix::t(A)
  seed.use <- abs(.Random.seed[[3]])
  set.seed(seed.use)
  if(is.matrix(k)){
    if(!(nrow(A) %in% dim(k))) stop("dimensions of matrix specified for 'k' are not compatible with number of rows in 'A'")
  } else {
    w_init <- matrix(runif(k * nrow(A)), k, nrow(A))
  }
  
  nmf_model <- c_gcnmf(A, At, G, tol, maxit, verbose, L1, L2, threads, w_init)
  rownames(nmf_model$h) <- colnames(nmf_model$w) <- paste0(reduction.key, 1:nrow(nmf_model$h))
  rownames(nmf_model$w) <- rnames
  colnames(nmf_model$h) <- cnames
  object@reductions[[reduction.name]] <- new("DimReduc",
                                             cell.embeddings = t(nmf_model$h),
                                             feature.loadings = nmf_model$w,
                                             assay.used = assay,
                                             stdev = nmf_model$d,
                                             global = FALSE,
                                             key = reduction.key)

  object
}

#' @rdname RunGCNMF
#'
#' @name RunGCNMF
#'
#' @export
#'
RunGCNMF <- function(object, ...) {
  UseMethod("RunGCNMF")
}

#' @rdname RunGCNMF
#'
#' @name RunGCNMF
#'
#' @export
#'
.S3method("RunGCNMF", "Seurat", RunGCNMF.Seurat)
