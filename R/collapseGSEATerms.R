#' Collapse redundant terms with set p-value
#'
#' @param object A Seurat, SingleCellExperiment, or RcppML::nmf object,
#'   a list of GSEA results, or a matrix of -log10 adjusted p-values.
#' @param reduction Dimensional reduction for which GSEA analysis has been performed.
#' @param gsea.name Name under which GSEA results are stored. Default is "gsea".
#' @param padj.sig Adjusted p-value threshold for selecting pathways to collapse.
#' @param pval.threshold p-value threshold for selecting significant pathways.
#' @param collapse.name name of collapsed pathways result
#' @param ... Additional arguments.
#'
#' @return object containing fgsea::collapsePathways result
#'
#' @export
collapseGSEATerms <- function(object, ...) {
  UseMethod("collapseGSEATerms")
}

#' @rdname collapseGSEATerms
#' @method collapseGSEATerms Seurat
#' @export
collapseGSEATerms.Seurat <- function(object,
                                     reduction = "nmf",
                                     gsea.name = "gsea",
                                     padj.sig = 0.05,
                                     pval.threshold = 0.05,
                                     collapse.name = "collapse",
                                     ...) {
  if (is.null(object@reductions[[reduction]]@misc[[gsea.name]])) {
    stop("GSEA results not found. Run RunGSEA() first.")
  }
  
  if (collapse.name %in% c("gsea", "ranks", "pathways")) {
    stop(
      "collapse.name cannot be one of the reserved names: ",
      "gsea, ranks, pathways.",
      call. = FALSE
    )
  }
  
  object@reductions[[reduction]]@misc[[gsea.name]][[collapse.name]] <-
    collapseGSEATerms(
      object = object@reductions[[reduction]]@misc[[gsea.name]], 
      padj.sig = padj.sig,
      pval.threshold = pval.threshold,
      ...
    )
  
  return(object)
}
.S3method("collapseGSEATerms", "Seurat", collapseGSEATerms.Seurat)

#' @rdname collapseGSEATerms
#' @method collapseGSEATerms SingleCellExperiment
#' @export
collapseGSEATerms.SingleCellExperiment <- function(object,
                                                   reduction = "nmf",
                                                   gsea.name = "gsea",
                                                   padj.sig = 0.05,
                                                   pval.threshold = 0.05,
                                                   collapse.name = "collapse",
                                                   ...) {
  
  if (collapse.name %in% c("gsea", "ranks", "pathways")) {
    stop(
      "collapse.name cannot be one of the reserved names: ",
      "gsea, ranks, pathways.",
      call. = FALSE
    )
  }
  
  object@metadata[[reduction]] <-
    collapseGSEATerms(
      object = object@metadata[[reduction]],
      gsea.name = gsea.name,
      padj.sig = padj.sig,
      pval.threshold = pval.threshold,
      collapse.name = collapse.name,
      ...
    )
  
  return(object)
}

.S3method(
  "collapseGSEATerms",
  "SingleCellExperiment",
  collapseGSEATerms.SingleCellExperiment
)

#' @rdname collapseGSEATerms
#' @method collapseGSEATerms nmfList
#' @export
collapseGSEATerms.nmfList <- function(object,
                                      gsea.name = "gsea",
                                      padj.sig = 0.05,
                                      pval.threshold = 0.05,
                                      collapse.name = "collapse",
                                      ...) {
  if (is.null(object[[gsea.name]])) {
    stop("GSEA results not found. Run RunGSEA() first.")
  }
  object[[gsea.name]][[collapse.name]] <-
    collapseGSEATerms(
      object = object[[gsea.name]],
      padj.sig = padj.sig,
      pval.threshold = pval.threshold, 
      ...
    )
  return(object)
}

.S3method("collapseGSEATerms", "nmfList", collapseGSEATerms.nmfList)

#' @rdname collapseGSEATerms
#' @method collapseGSEATerms nmf
#' @export
collapseGSEATerms.nmf <- function(object,
                                  gsea.name = "gsea",
                                  padj.sig = 0.05,
                                  pval.threshold = 0.05,
                                  collapse.name = "collapse",
                                  ...) {
  if (is.null(object@misc[[gsea.name]])) {
    stop("GSEA results not found. Run RunGSEA() first.")
  }
  object@misc[[gsea.name]][[collapse.name]] <-
    collapseGSEATerms(
      object = object@misc[[gsea.name]],
      padj.sig = padj.sig,
      pval.threshold = pval.threshold,
      ...
    )
  
  
  return(object)
}

.S3method("collapseGSEATerms", "nmf", collapseGSEATerms.nmf)

#' @rdname collapseGSEATerms
#' @method collapseGSEATerms list
#' @importFrom fgsea collapsePathways
#' @export
collapseGSEATerms.list <- function(object,
                                   padj.sig = 0.05,
                                   pval.threshold = 0.05,
                                   ...) {
  if (!all(c("gsea", "ranks", "pathways") %in% names(object))) {
    stop("GSEA results not found. Run RunGSEA() first.")
  }
  
  collapsed =  collapseGSEA(
    gsea = object, 
    padj.sig = padj.sig,
    pval.threshold = pval.threshold,
    ...)
  
  
  return(collapsed)
}
.S3method("collapseGSEATerms", "list", collapseGSEATerms.list)
