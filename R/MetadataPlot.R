#' @rdname RunLNMF
#' @name MetadataPlot
#' @export
#'
MetadataPlot.Seurat <- function(object, split.by, reduction = "lnmf", ...) {
  if (!(reduction %in% names(object@reductions))) {
    stop("this Seurat object does not contain the requested reductions slot")
  }
  plot(MetadataSummary(t(object@reductions[[reduction]]@cell.embeddings), object@meta.data[[split.by]]))
}


#' @rdname RunLNMF
#' @export
#'
.S3method("MetadataPlot", "Seurat", MetadataPlot.Seurat)


#' @rdname RunLNMF
#' @export
#'
MetadataPlot <- function(object, ...) {
  UseMethod("MetadataPlot")
}
