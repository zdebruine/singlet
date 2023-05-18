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
  if (is.null(object@assays[[assay]]@key)) {
    object@assays[[assay]]@key <- paste0(assay, "_")
  }
  object@assays[[assay]] <- PreprocessData(object@assays[[assay]], ...)
  object
}

#' @rdname PreprocessData
#' @export
PreprocessData.Assay <- function(object, scale.factor = 10000, ...) {
  if (ncol(object@counts) == 0) {
    object@data <- PreprocessData(object@data, ...)
  } else {
    object@data <- PreprocessData(object@counts, ...)
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
