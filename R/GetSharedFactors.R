#' @export
#' @rdname RunLNMF
#'
GetSharedFactors <- function(object, split.by, reduction = "lnmf") {
  if (!(reduction %in% names(object@reductions))) {
    stop("this Seurat object does not contain the requested reductions slot")
  }
  # which(rowSums(object@reductions[[reduction]]@misc$link_matrix == 0) == 0)
  which(!(colnames(object@reductions[[reduction]]@cell.embeddings) %in% names(which(apply(MetadataSummary(t(object@reductions[[reduction]]@cell.embeddings), object@meta.data[[split.by]]), 2, function(x) min(x) == 0)))))
}
