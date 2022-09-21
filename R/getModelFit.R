#' get linear all-pairs comparisons fits for a design matrix and data matrix
#' 
#' Continuing along with the theme of "stupid limma tricks", this function 
#' fits and (optionally, though by default) shrinks a linear model for a factor.
#'
#' @param design      a model.matrix (or a sparse.model.matrix, perhaps)
#' @param object      a data.matrix or a Seurat DimReduc object
#' @param shrink      run limma::eBayes() on the model? (TRUE) 
#'
#' @export
getModelFit <- function(design, object, shrink=TRUE) {

  dat <- object
  if (is(object, "nmf")) dat <- object@h
  if (is(object, "DimReduc")) dat <- t(object@cell.embeddings) # DimReduc
  stopifnot(all(rownames(design) %in% colnames(dat)))
  fit <- lmFit(dat[, rownames(design)], design)
  if (shrink) fit <- eBayes(fit) 
  return(fit)

}
