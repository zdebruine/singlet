#' Project a factor model
#'
#' @description Project a dataset onto a factor model for transfer learning
#'
#' @inheritParams run_nmf
#' @param w matrix giving the factor model, of dimensions \code{nrow(A) x k}
#' @return list of \code{h} and \code{d}, where \code{d} gives the relative contribution of each factor in \code{h} to the model
#'
#' @export
#'
project_model <- function(A, w, L1 = 0.01, L2 = 0, threads = 0) {
  if (nrow(w) != nrow(A) & 
      ncol(w) != nrow(A)) {
    stop("'w' must share a common edge with the rows of 'A'")
  } else { 
    # leftover from testing 
    singlet:::c_project_model(A, w, L1, L2, threads)
  }
}


#' Project data onto a factor model
#'
#' @description Non-negative Least Squares (NNLS) projection of assay data onto a factor model for transfer learning
#'
#' @inheritParams RunNMF.Seurat
#' @param w a factor model of the same number of rows as the assay to be projected in the current object, where columns correspond to factors
#' @param reorder reorder the factors of the projection by d? (FALSE) 
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
                          reorder = FALSE, 
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

  if (reorder) { 
    idx <- order(nmf_model$d, decreasing = TRUE)
    nmf_model$h <- nmf_model$h[idx, ]
    nmf_model$d <- nmf_model$d[idx]
    nmf_model$w <- nmf_model$w[, idx]
  }
  
  object@reductions[[reduction.name]] <- new("DimReduc",
                                             cell.embeddings = t(nmf_model$h),
                                             feature.loadings.projected = nmf_model$w,
                                             assay.used = assay,
                                             stdev = nmf_model$d,
                                             key = reduction.key)
  object
}


#' Project data onto a factor model
#'
#' @description Non-negative Least Squares (NNLS) projection of assay data onto a factor model for transfer learning
#'
#' @inheritParams RunNMF.SingleCellExperiment
#'
#' @param w       factor loadings with nrow(w) equal to nrow(object)
#' @param reorder reorder the factors of the projection by d? (FALSE) 
#' @param ...     not implemented
#'
#' @details Use \code{set.seed()} to guarantee reproducibility!
#'
#' @aliases ProjectData
#'
#' @seealso \code{\link{RunLNMF}}, \code{\link{MetadataSummary}}
#'
#' @rdname ProjectData
#'
#' @return a SingleCellExperiment with projection stored in reducedDim(, "NNLS")
#'
#' @import scuttle
#'
#' @export
ProjectData.SingleCellExperiment <- function(object, 
                                             w,
                                             split.by = NULL,
                                             assay = "logcounts",
                                             L1 = 0.01,
                                             L2 = 0,
                                             reduction.name = "NMF",
                                             reduction.key = "NMF_",
                                             threads = 0,
                                             reorder = FALSE,
                                             ...) {

  # check if data exists and has been normalized
  if (assay == "logcounts" && !assay %in% assayNames(object)) {
    object <- logNormCounts(object)
  }

  # pull the data for projection 
  A <- as(assay(object, assay), "CsparseMatrix")
  rnames <- rownames(object)
  cnames <- colnames(object)
 
  # reweight? 
  if (!is.null(split.by)) {
    split.by <- as.integer(as.numeric(as.factor(colData(object)[,split.by])))-1
    if (any(sapply(split.by, is.na))) stop("'split.by' cannot contain NAs")
    A <- weight_by_split(A, split.by, length(unique(split.by)))
  }
  
  # check if we can proceed 
  stopifnot(assay %in% assayNames(object))
  if (!all(rownames(x) %in% rownames(A)) | 
      !all(rownames(A) %in% rownames(w))) {
    warning("Some genes are missing from the data to be projected.")
  }
  w <- w[which(rownames(w) %in% rownames(A)), ]
  A <- A[which(rownames(A) %in% rownames(w)), ]
 
  # project
  nmf_model <- project_model(A, w, L1, L2, threads)
  nmf_model$w <- w

  # label
  colnames(nmf_model$w) <- paste0(reduction.key, 1:nrow(nmf_model$h))
  rownames(nmf_model$h) <- colnames(nmf_model$w)
  rownames(nmf_model$w) <- rownames(A)
  colnames(nmf_model$h) <- cnames
 
  # reorder?
  if (reorder) { 
    idx <- order(nmf_model$d, decreasing = TRUE)
    nmf_model$h <- nmf_model$h[idx, ]
    nmf_model$d <- nmf_model$d[idx]
    nmf_model$w <- nmf_model$w[, idx]
  }
  
  # store the results
  lem <- LinearEmbeddingMatrix(sampleFactors = t(nmf_model$h),
                               featureLoadings = nmf_model$w,
                               metadata = list(d = nmf_model$d))
  reducedDim(object, reduction.name) <- lem

  # return
  object

}


#' @export
#' @rdname ProjectData
#'
ProjectData <- function(object, ...) {
  UseMethod("ProjectData")
}
