#' verify that a list of matrices is in fact a named list of model matrices
#' 
#' @param designs     an alleged list of model matrices
#' 
#' @return            the list of model matrices, assuming it passes
#'
#' @details           this function will squawk and stop if the list is no good
#'
#' @export
checkDesigns <- function(designs) {

  if (is.null(names(designs)) |
      !all(sapply(designs, function(x) !is.null(attr(x, "assign"))))) {
    stop("`designs` must be a named list of model.matrix outputs.")
  } else { 
    return(designs)
  }

}
