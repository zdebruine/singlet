#' automatically generate a means model (one-vs-all group associations) 
#'
#' A little-known trick in limma is to fit ~ 0 + group for a means model.
#' This function automates that for a data.frame and a factor column of it.
#'
#' @param field       the name of a column in the data.frame, or the column
#' @param meta.data   a data.frame with one or more factor columns, or NULL
#' @param sparse      fit a sparse model.matrix? (FALSE) 
#' @param ova         fit a One-Vs-All model matrix (no referent)? (TRUE)
#' @param ...         any additional params to pass to model.matrix
#'
#' @details
#' If a factor (and no meta.data) is supplied (usually by with(meta.data, ...)),
#' getModelMatrix will attempt to figure out the text to remove from the matrix
#' column names by using deparse() and match.call() on the arguments (voodoo!).
#' In order to fit one-vs-all comparisons, a means model is the default. If you
#' have a referent group (e.g. normal bone marrow vs. a bunch of leukemia cells)
#' or simply don't want a means model, set `ova` (one vs all) to FALSE. 
#' 
#' @return            a model.matrix or sparse.model.matrix (if sparse==TRUE)
#'
#' @import Matrix
#'
#' @export
getModelMatrix <- function(field, meta.data=NULL, sparse=FALSE, ova=TRUE, ...) {

  if (is.null(meta.data)) { 
    if (is.factor(field) & nlevels(field) > 1) {
      fieldname <- as.character(sapply(match.call()[-1], deparse)[1]) # voodoo
      meta.data <- data.frame(field)
      names(meta.data) <- fieldname
      field <- fieldname
    } else { 
      stop("If meta.data is NULL, `field` must be a factor with > 1 levels.")
    }
  } else {
    stopifnot(field %in% names(meta.data))
  }

  notNA <- which(is.na(meta.data[[field]]))
  if (!sparse) {
    if (ova) {
      mat <- model.matrix(~ 0 + ., data=meta.data[, field, drop=FALSE], ...)
    } else { 
      message("Fitting a model with a referent group. Be sure you want this!")
      mat <- model.matrix(~ ., data=meta.data[, field, drop=FALSE], ...)
    }
  } else { 
    if (ova) { 
      mat <- sparse.model.matrix(~ 0 + ., data=meta.data[, field, drop=FALSE])
    } else {
      message("Fitting a model with a referent group. Be sure you want this!")
      mat <- sparse.model.matrix(~ ., data=meta.data[, field, drop=FALSE], ...)
    }
  }
  colnames(mat) <- gsub(field, "", colnames(mat))
  return(mat)

}
