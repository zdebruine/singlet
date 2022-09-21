#' Project a factor model
#' 
#' @description Project a dataset onto a factor model for transfer learning
#' 
#' @inheritParams run_nmf
#' @param w matrix giving the factor model, of dimensions \code{nrow(A) x k}
#' @return list of \code{h} and \code{d}, where \code{d} gives the relative contribution of each factor in \code{h} to the model
#' @export
project_model <- function(A, w, L1 = 0.01, L2 = 0, threads = 0){
  if(nrow(w) != nrow(A) & ncol(w) != nrow(A)) stop("'w' must share a common edge with the rows of 'A'")

  c_project_model(A, w, L1, L2, threads)
}

