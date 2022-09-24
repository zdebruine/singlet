#' annotate an NMF model 
#' 
#' @param object    an object suitable for annotation (Seurat, DimReduc, or nmf)
#' @param fields    fields (columns) of meta.data to annotate against
#' @param meta.data a data.frame
#'
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
#' @inheritParams AnnotateNMF
#' @param shrink    apply eBayes to the fits? (TRUE) 
#'
#' @rdname AnnotateNMF
#' @aliases AnnotateNMF
#'
#' @import limma
#'
#' @examples 
#' if (!exists("pbmc3k")) get_pbmc3k_data() %>% NormalizeData -> pbmc3k
#' if (!"nmf" %in% Reductions(pbmc3k)) pbmc3k %>% RunNMF() -> pbmc3k
#' AnnotateNMF(pbmc3k)
#' 
#' @export
#'
AnnotateNMF.DimReduc <- function(object, meta.data, shrink=TRUE, ...){

  # get all factors with multiple levels from meta.data data.frame
  eligible <- sapply(meta.data, function(x) is(x, "factor") & nlevels(x) > 1)
  if (sum(eligible) < 1) stop("No factors in meta.data slot with >1 level!")
  fields <- names(which(eligible))
  names(fields) <- fields

  designs <- lapply(fields, getModelMatrix, meta.data=meta.data)
  fits <- lapply(designs, getModelFit, object=object, shrink=shrink)
  annotations <- lapply(fits, getModelResults)

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
#'
#' @inheritParams AnnotateNMF
#'
#' @param reduction the reductions slot in the Seurat object containing the model to annotate
#'
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


#' Annotate NMF model with cell metadata
#' 
#' Note: some of the documentation (e.g. meta.data is a named list) is wrong.
#' 
#' @details Maps factor information in an RcppML::nmf object against meta.data
#' 
#' @inheritParams AnnotateNMF
#' @param shrink    apply eBayes to the fits? (TRUE) 
#'
#' @rdname AnnotateNMF
#' @aliases AnnotateNMF
#'
#' @import limma
#'
#' @export
#'
AnnotateNMF.nmf <- function(object, meta.data, shrink=TRUE, ...){

  # get all factors with multiple levels from meta.data data.frame
  eligible <- sapply(meta.data, function(x) is(x, "factor") & nlevels(x) > 1)
  if (sum(eligible) < 1) stop("No factors in meta.data slot with >1 level!")
  fields <- names(which(eligible))
  names(fields) <- fields

  designs <- lapply(fields, getModelMatrix, meta.data=meta.data)
  # it turns out that we can, and did, refactor Seurat vs RcppML::nmf into a fn
  fits <- lapply(designs, getModelFit, object=object, shrink=shrink)
  annotations <- lapply(fits, getModelResults)

  # we use this slot in RcppML::nmf same as in Seurat DimReduc
  object@misc$annotations <- annotations
  object

}


#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
#' @export
#'
.S3method("AnnotateNMF", "nmf", AnnotateNMF.nmf)
