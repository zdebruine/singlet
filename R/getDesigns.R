#' Refactored out from AnnotateNMF to ease argument handling
#'
#' @param columns   factor columns of meta.data, optional if !is.null(designs)
#' @param meta.data a data.frame of annotations, optional if !is.null(designs)
#' @param designs   named list of design matrices (supersedes meta.data/columns)
#' @param max.levels maximum number of levels permitted for a factor to be kept
#'
#' @return a named list of design matrices, if one was not provided
#' @export
getDesigns <- function(columns = NULL, meta.data = NULL, designs = NULL, max.levels = 200) {

  if (is.null(designs)) {
    stopifnot(any(!is.null(c(columns, meta.data))))
    columns <- checkColumns(meta.data = meta.data, 
                            columns = columns, 
                            max.levels = max.levels)
    designs <- lapply(columns, getModelMatrix, meta.data = meta.data)
  }

  checkDesigns(designs)
}
