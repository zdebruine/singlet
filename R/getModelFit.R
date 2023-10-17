#' get linear all-pairs comparisons fits for a design matrix and data matrix
#' 
#' Continuing along with the theme of "stupid limma tricks", this function 
#' fits and shrinks a means model for a factor. The proportion of factors
#' assumed to have a fold-change > 0 is 1%, and a robust fit is applied. 
#'
#' @param design      a model.matrix (or a sparse.model.matrix, perhaps)
#' @param object      a data.matrix, Seurat DimReduc, or RcppML nmf object
#' @param center      center the factor matrix for testing? (TRUE) 
#' @param ...         additional arguments, passed to base::scale
#'
#' 
#' @examples 
#' if (FALSE) { 
#'   get_pbmc3k_data() %>% NormalizeData() -> pbmc3k
#'   design <- model.matrix(~ 0 + cell_type, data=pbmc3k@meta.data)
#'   fit <- getModelFit(design, pbmc3k) # toy fit on lognormcounts  
#'   # Subsetting data to non-NA observations to match design matrix.
#'   limma::topTable(fit)
#' }
#' 
#' @export
getModelFit <- function(design, object, center=TRUE, ...) {

  dat <- object
  if (is(object, "nmf")) dat <- object@h # RcppML nmf 
  if (is(object, "Seurat")) dat <- object@assays$RNA@data
  if (is(object, "DimReduc")) dat <- t(object@cell.embeddings)
  if (is(object, "SingleCellExperiment")) dat <- logcounts(object)
  # SingleCellExperiment::reducedDim(object, dimname) just returns a data.matrix
  if (nrow(dat) < nrow(design)) dat <- t(dat) # transpose reduced dims if needed

  # janky, but should be foolproof 
  if (nrow(design) != nrow(dat)) {
    if (!all(rownames(design) %in% colnames(dat))) {
      message("Rows of the design matrix do not match columns of the object.")
      message("This usually means that there are NAs in the sample metadata.")
      message("Ensure rownames of your design matrix match data column names.")
      message("Alternatively, provide object[, !is.na(object$predictor)]")
      message("so that the dimensions of the data and design matrices match.")
      stop("Cannot proceed as called.")
    } else { 
      message("Subsetting data to non-NA observations to match design matrix.")
      tofit <- dat[, rownames(design)]
    }
  } else {
    tofit <- dat
    if (is.null(rownames(design))) {
      warning("Design matrix has appropriate rank, but no row names. Beware!")
    } else if (!identical(rownames(design), colnames(tofit))) {
      warning("Design matrix row names do not match data observation names!")
      warning("This is usually a VERY BAD THING. You MUST check your data.")
      warning("If this warning message persists, file a bug with a reprex.")
    }
  }

  if (center) tofit <- t(scale(t(tofit), ...))
  fit <- eBayes(lmFit(tofit, design), proportion=0.01, robust=TRUE)
  fit$centered <- center
  return(fit)

}
