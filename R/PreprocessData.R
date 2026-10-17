#' Normalize count data
#'
#' Standard log-normalization equivalent to \code{Seurat::LogNormalize}
#'
#' @param object Seurat object
#' @param assay assay in which the counts matrix resides
#' @param scale.factor value by which to multiply all columns after unit normalization and before \code{log1p} transformation
#' @param ... arguments to \code{Seurat::LogNormalize}
#' @export
#' @rdname PreprocessData
#'
PreprocessData.Seurat <- function(object, scale.factor = 10000, assay = NULL, ...) {
  if (is.null(assay)) assay <- names(object@assays)[[1]]
  # Normalize the raw counts when the assay has them, otherwise whatever is in the data layer.
  # This goes through the layer accessors, so it serves Assay (v3/v4) and Assay5 (v5) alike.
  counts <- .assay_matrix(object, assay, "counts")
  if (is.null(counts)) counts <- .assay_matrix(object, assay, "data")
  if (is.null(counts)) {
    stop("assay '", assay, "' has neither a 'counts' nor a 'data' layer to normalize. ",
         "If the layers are split, join them first with SeuratObject::JoinLayers().", call. = FALSE)
  }
  .set_assay_matrix(object, assay, "data", PreprocessData(counts, scale.factor = scale.factor, ...))
}

#' @rdname PreprocessData
#' @export
PreprocessData.Assay <- function(object, scale.factor = 10000, ...) {
  if (ncol(object@counts) == 0) {
    object@data <- PreprocessData(object@data, scale.factor = scale.factor, ...)
  } else {
    object@data <- PreprocessData(object@counts, scale.factor = scale.factor, ...)
  }
  object
}

#' @rdname PreprocessData
#' @export
PreprocessData.dgCMatrix <- function(object, scale.factor = 10000, ...) {
  m <- Seurat::LogNormalize(object, scale.factor, ...)
  rownames(m) <- rownames(object)
  colnames(m) <- colnames(object)
  m
}

#' @export
#' @rdname PreprocessData
#'
PreprocessData <- function(object, scale.factor, ...) {
  UseMethod("PreprocessData")
}

#' @export
#' @rdname PreprocessData
#' @name PreprocessData
#'
.S3method("PreprocessData", "dgCMatrix", PreprocessData.dgCMatrix)

#' @export
#' @rdname PreprocessData
#' @name PreprocessData
#'
.S3method("PreprocessData", "Assay", PreprocessData.Assay)


#' @export
#' @rdname PreprocessData
#' @name PreprocessData
#'
.S3method("PreprocessData", "Seurat", PreprocessData.Seurat)
