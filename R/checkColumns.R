#' verify that columns for auto-annotation are factors with > 1 level
#'
#' @param meta.data the meta.data (or a Seurat object if needs be)
#' @param columns   the columns (optional; if NULL, will check all columns)
#' @param max.levels maximum number of levels permitted for a factor to be kept
#' @return          a vector of suitable columns (may be length 0)
#'
#' @export
checkColumns <- function(meta.data, columns = NULL, max.levels = 200) {
  verbose <- !is.null(columns)
  if (is(meta.data, "Seurat")) meta.data <- meta.data@meta.data
  if (is.null(columns)) columns <- colnames(meta.data)
  names(columns) <- columns
  keep <- names(which(sapply(columns, .keepColumn, meta.data = meta.data, max.levels = max.levels)))
  discard <- setdiff(columns, keep)
  if (verbose & length(discard) > 0) {
    message("Some columns are not factors, or have only one level, or have more than max.levels levels.")
    message("Skipping `", paste(discard, collapse = "`, `"), "`.")
  }
  names(keep) <- keep
  return(keep)
}


# helper fn
.keepColumn <- function(x, meta.data, max.levels) {
  if (!x %in% names(meta.data)) {
    return(FALSE)
  }
  if (!is(meta.data[[x]], "factor")) {
    return(FALSE)
  }
  if (nlevels(meta.data[[x]]) < 2) {
    return(FALSE)
  }
  if (nlevels(meta.data[[x]]) > max.levels) {
    return(FALSE)
  }
  return(TRUE)
}
