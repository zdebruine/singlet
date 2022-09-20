#' @export
#'
AnnotateNMF <- function(object, ...) {
  UseMethod("AnnotateNMF")
}


#' Annotate NMF model with cell metadata
#' 
#' Note: some of the documentation (e.g. meta.data is a named list) is wrong.
#' This function clearly expects meta.data to be a data.frame with cells as rows
#' 
#' @details Maps factor information from the \code{meta.data} slot of a Seurat object or a specified factor to each NMF factor, and computes summary statistics
#' 
#' @inheritParams RunNMF
#' @param fields if a Seurat object is specified, fields in the \code{meta.data} slot may be specified for which to compute summary statistics. If \code{fields = NULL}, all fields that are a \code{factor} with more than one level are considered.
#' @param meta.data if a DimReduc object is specified, provide a named list containing factor(s) for which to calculate summary statistics against the NMF model
#'
#' @param shrink  apply eBayes to the fits? (TRUE) 
#'
#' @rdname AnnotateNMF
#' @aliases AnnotateNMF
#'
#' @import limma
#'
#' @export
#'
AnnotateNMF.DimReduc <- function(object, meta.data, shrink=TRUE, ...){

  # get all factors with multiple levels from meta.data data.frame
  eligible <- sapply(meta.data, function(x) is(x, "factor") & nlevels(x) > 1)
  if (sum(eligible) < 1) stop("No factors in meta.data slot with >1 level!")
  fields <- names(which(eligible))
  names(fields) <- fields


  # create means models for each field (can also use sparse.model.matrix)
  getModelMatrix <- function(field, meta.data) {
    notNA <- which(is.na(meta.data[, field]))
    mat <- model.matrix(~ 0 + ., data=meta.data[, field, drop=FALSE])
    colnames(mat) <- gsub(field, "", colnames(mat))
    return(mat)
  }
  designs <- lapply(fields, getModelMatrix, meta.data=meta.data)


  # get linear all-pairs comparisons fits for each field
  getModelFit <- function(design, object, shrink=TRUE) {
    fit <- lmFit(t(object@cell.embeddings[rownames(design), ]), design)
    if (shrink) fit <- eBayes(fit) 
    return(fit)
  }
  fits <- lapply(designs, getModelFit, object=object, shrink=shrink)


  # log-odds for differential representation of a factor in each group can be 
  # found in fit$lods (for any given fit, this will be a matrix)
  getGroupResults <- function(fit) { 

    lods <- reshape2::melt(fit$lods)
    names(lods) <- c("factor", "group", "fc")
    pBH <- reshape2::melt(apply(fit$p.value, 2, p.adjust, method="BH"))
    names(pBH) <- c("factor", "group", "p")
    merge(lods, pBH)

  } 

  annotations <- lapply(fits, getGroupResults)
  object@misc$annotations <- annotations
  object
}


#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
#' @export
#'
.S3method("AnnotateNMF", "DimReduc", AnnotateNMF.DimReduc)


#' @rdname AnnotateNMF
#' @param reduction the reductions slot in the Seurat object containing the model to annotate
#' @aliases AnnotateNMF
#'
#' @export
#' 
AnnotateNMF.Seurat <- function(object, fields = NULL, reduction = "nmf", ...){

  if(!is.null(fields)){
    for(i in 1:length(fields)){
      if(!is.factor(object@meta.data[[fields[[i]]]])){
        stop(fields[[i]], "is not a factor!")
      } else if(length(levels(object@meta.data[[fields[[i]]]]))){
        stop(fields[[i]], "is a factor, but only has one level")
      }
    }

    object@reductions[[reduction]] <- AnnotateNMF.DimReduc(object@reductions[[reduction]], object@meta.data[, fields], ...)

  } else {

    object@reductions[[reduction]] <- AnnotateNMF.DimReduc(object@reductions[[reduction]], object@meta.data)

  }

  object

}


#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
#' @export
#'
.S3method("AnnotateNMF", "Seurat", AnnotateNMF.Seurat)
