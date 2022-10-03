#' get linear all-pairs comparisons fits for a design matrix and data matrix
#' 
#' Continuing along with the theme of "stupid limma tricks", this function 
#' fits and shrinks a means model for a factor. The proportion of factors
#' assumed to have a fold-change > 0 is 1%, and a robust fit is applied. 
#'
#' @param design      a model.matrix (or a sparse.model.matrix, perhaps)
#' @param object      a data.matrix, Seurat DimReduc, or RcppML nmf object
#'
#' @export
getModelFit <- function(design, object) {

  dat <- object
  if (is(object, "nmf")) dat <- object@h # RcppML nmf 
  if (is(object, "DimReduc")) dat <- t(object@cell.embeddings) # Seurat DimReduc
  # SingleCellExperiment::reducedDim(object, dimname) just returns a data.matrix
  stopifnot(all(rownames(design) %in% colnames(dat)))
  dat <- dat[, rownames(design)]

  # center the matrix,
  # to force exclusivity
  rmeans <- rowMeans(dat)
  names(rmeans) <- rownames(dat)
  tofit <- sweep(dat, 1, rmeans, `-`)
  fit <- eBayes(lmFit(tofit, design), proportion=0.01, robust=TRUE)
  fit$rowmeans <- rmeans
  return(fit)

}
