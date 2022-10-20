#' verify that columns for auto-annotation are factors with > 1 level
#' 
#' @param meta.data the meta.data (or a Seurat object if needs be)
#' @param columns   the columns (optional; if NULL, will check all columns)
#'
#' @return          a vector of suitable columns (may be length 0)
#'
#' @export
checkColumns <- function(meta.data, columns = NULL) {

  verbose <- !is.null(columns)
  if (is(meta.data, "Seurat")) meta.data <- meta.data@meta.data
  if (is.null(columns)) columns <- colnames(meta.data)
  names(columns) <- columns
  keep <- names(which(sapply(columns, .keepColumn, meta.data = meta.data)))
  discard <- setdiff(columns, keep)
  if (verbose & length(discard) > 0) {
    warning("Some columns are not factors, or have only one level.")
    warning(paste("Skipping ",  paste(discard, collapse=", ")))
  }
  names(keep) <- keep
  return(keep)

}


# helper fn
.keepColumn <- function(x, meta.data) {
  
  if (!x %in% names(meta.data)) return(FALSE)
  if (!is(meta.data[[x]], "factor")) return(FALSE)
  if (nlevels(meta.data[[x]]) < 2) return(FALSE) 

  return(TRUE) 

}
