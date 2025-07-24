#' annotate an NMF model
#'
#' @param object    an object suitable for annotation (Seurat, DimReduc, or nmf)
#' @param columns   factor columns of meta.data (see below) to annotate against
#' @param meta.data a data.frame, if one is not already part of the object
#' @param designs   named list of design matrices (supersedes meta.data/columns)
#' @param center    center the factor matrix for testing? (TRUE)
#' @param scale     scale the factor matrix for testing? (FALSE)
#' @param max.levels maximum number of levels a factor may have in order to be included in analysis
#' @param tail      add direction for one or two-tailed testing (pos,neg,std), default = pos
#' @param annotation.name name of factor annotation, default = "annotations"
#' @param ... not implemented
#' @export
#'
AnnotateNMF <- function(object, ...) {
  UseMethod("AnnotateNMF")
}


#' Annotate NMF model with cell or sample metadata
#'
#' @rdname AnnotateNMF
#' @aliases AnnotateNMF
#'
#' @import limma
#'
#' @export
#'
AnnotateNMF.DimReduc <- function(object, meta.data = NULL, columns = NULL, designs = NULL, center = TRUE, scale = FALSE, max.levels = 200, tail = "pos", annotation.name = "annotations",...) {
  designs <- getDesigns(columns = columns,
                        meta.data = meta.data, 
                        designs = designs, 
                        max.levels)
  fits <- lapply(designs, 
                 getModelFit, 
                 object = object, 
                 center = center, 
                 scale = scale)
  object@misc[[annotation.name]] <- lapply(fits, getModelResults, tail = tail)
  return(object)
}


#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
#' @export
#'
.S3method("AnnotateNMF", "DimReduc", AnnotateNMF.DimReduc)


#' @rdname AnnotateNMF
#'
#' @param reduction the reductions slot in the Seurat object containing the model to annotate
#'
#' @examples
#' \dontrun{
#' get_pbmc3k_data() %>%
#'   NormalizeData() %>%
#'   RunNMF() -> pbmc3k
#' AnnotateNMF(pbmc3k)
#' }
#' @aliases AnnotateNMF
#'
#' @export
#'
AnnotateNMF.Seurat <- function(object, columns = NULL, reduction = "nmf", tail = "pos", annotation.name = "annotations",...) {
  if (is.null(columns)) columns <- colnames(object@meta.data)
  object@reductions[[reduction]] <-
    AnnotateNMF.DimReduc(
      object = object@reductions[[reduction]],
      meta.data = object@meta.data[, columns],
      columns = columns, 
      tail = tail,
      annotation.name = annotation.name,...
    )
  return(object)
}


#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
#' @export
#'
.S3method("AnnotateNMF", "Seurat", AnnotateNMF.Seurat)


#' Annotate NMF model with cell metadata
#'
#' @details Maps factor information in an RcppML::nmf object against meta.data
#'
#' @rdname AnnotateNMF
#' @aliases AnnotateNMF
#'
#' @import limma
#'
#' @export
#'
AnnotateNMF.nmf <- function(object, meta.data, columns = NULL, designs = NULL, center = TRUE, scale = FALSE, max.levels = 200, tail = "pos",annotation.name = "annotations", ...) {
  designs <- getDesigns(columns = columns, meta.data = meta.data, designs = designs, max.levels, ...)
  fits <- lapply(designs, getModelFit, object = object, center = center, scale = scale)
  object@misc[[annotation.name]] <- lapply(fits, getModelResults, tail = tail)
  return(object)
}


#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
#' @export
#'
.S3method("AnnotateNMF", "nmf", AnnotateNMF.nmf)
