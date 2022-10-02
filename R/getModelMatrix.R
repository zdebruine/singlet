#' automatically generate a means model (one-vs-all group associations) 
#'
#' A little-known trick in limma is to fit ~ 0 + group for a means model.
#' This function automates that for a data.frame and a factor column of it.
#'
#' @param field       the name of a column in the data.frame
#' @param meta.data   a data.frame with one or more factor columns
#' @param sparse      fit a sparse model.matrix? (FALSE) 
#' @param ...         any additional params to pass to model.matrix
#'
#' @return            a model.matrix or sparse.model.matrix (if sparse==TRUE)
#'
#' @examples 
#' 
#' if (!exists("pbmc3k") | ! "nmf" %in% Reductions(pbmc3k)) { 
#'   get_pbmc3k_data() %>% NormalizeData %>% RunNMF %>% AnnotateNMF -> pbmc3k
#' } 
#' 
#'
#' @import Matrix
#'
#' @export
getModelMatrix <- function(field, meta.data, sparse=FALSE, ...) {

  stopifnot(field %in% names(meta.data))
  notNA <- which(is.na(meta.data[, field]))
  if (!sparse) {
    mat <- model.matrix(~ 0 + ., data=meta.data[, field, drop=FALSE], ...)
  } else { 
    mat <- sparse.model.matrix(~ 0 + ., data=meta.data[, field, drop=FALSE])
  }
  colnames(mat) <- gsub(field, "", colnames(mat))
  return(mat)

}
