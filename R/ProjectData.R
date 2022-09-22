#' Project a factor model
#'
#' @description Project a dataset onto a factor model for transfer learning
#'
#' @inheritParams run_nmf
#' @param w matrix giving the factor model, of dimensions \code{nrow(A) x k}
#' @return list of \code{h} and \code{d}, where \code{d} gives the relative contribution of each factor in \code{h} to the model
#' @export
project_model <- function(A, w, L1 = 0.01, L2 = 0, threads = 0) {
  if (nrow(w) != nrow(A) & ncol(w) != nrow(A)) stop("'w' must share a common edge with the rows of 'A'")

  c_project_model(A, w, L1, L2, threads)
}


#' Project data onto a factor model
#'
#' @description Non-negative Least Squares (NNLS) projection of assay data onto a factor model for transfer learning
#'
#' @inheritParams RunNMF.Seurat
#' @param w a factor model of the same number of rows as the assay to be projected in the current object, where columns correspond to factors
#' @param ... not implemented
#' @details Use \code{set.seed()} to guarantee reproducibility!
#' @export
#' @aliases ProjectData
#' @seealso \code{\link{RunLNMF}}, \code{\link{MetadataSummary}}
#' @rdname ProjectData
#' @return Returns a Seurat object with the projection stored in the reductions slot
#'
ProjectData.Seurat <- function(object, w,
                          split.by = NULL,
                          assay = NULL,
                          L1 = 0.01,
                          L2 = 0,
                          reduction.name = "nmf_projection",
                          reduction.key = "NNLS_",
                          threads = 0,
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
  rnames <- rownames(A)
  cnames <- colnames(A)
  
  if (!is.null(split.by)) {
    split.by <- as.integer(as.numeric(as.factor(object@meta.data[[split.by]]))) - 1
    if (any(sapply(split.by, is.na))) {
      stop("'split.by' cannot contain NA values")
    }
    A <- weight_by_split(A, split.by, length(unique(split.by)))
  }
  
  w <- w[which(rownames(w) %in% rownames(A)), ]
  A <- A[which(rownames(A) %in% rownames(w)), ]
  
  nmf_model <- project_model(A, w, L1, L2, threads)
  nmf_model$w <- w
  rownames(nmf_model$h) <- colnames(nmf_model$w) <- paste0(reduction.key, 1:nrow(nmf_model$h))
  rownames(nmf_model$w) <- rownames(A)
  colnames(nmf_model$h) <- cnames
  idx <- order(nmf_model$d, decreasing = TRUE)
  nmf_model$h <- nmf_model$h[idx, ]
  nmf_model$d <- nmf_model$d[idx]
  nmf_model$w <- nmf_model$w[, idx]
  object@reductions[[reduction.name]] <- new("DimReduc",
                                             cell.embeddings = t(nmf_model$h),
                                             feature.loadings.projected = nmf_model$w,
                                             assay.used = assay,
                                             stdev = nmf_model$d,
                                             key = reduction.key)
  object
}


#' @export
#' @rdname ProjectData
#'
ProjectData <- function(object, ...) {
  UseMethod("ProjectData")
}
