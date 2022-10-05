#' refactored out from AnnotateNMF to ease argument handling
#'
#' @param object    a suitable data structure for validating design matrices
#' @param columns   factor columns of meta.data, optional if !is.null(designs)
#' @param meta.data a data.frame of annotations, optional if !is.null(designs)
#' @param designs   named list of design matrices (supersedes meta.data/columns)
#'
#' @return          a named list of design matrices, if one was not provided
#'
#' @export
getDesigns <- function(columns=NULL, meta.data=NULL, designs=NULL) {

  if (is.null(designs)) {
    stopifnot(any(!is.null(c(columns, meta.data))))
    columns <- checkColumns(meta.data=meta.data, columns=columns) 
    designs <- lapply(columns, getModelMatrix, meta.data=meta.data)
  }

  checkDesigns(designs)

}
