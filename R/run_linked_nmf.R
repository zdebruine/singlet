#' Run Linked Non-negative Matrix Factorization
#'
#' @description Run LNMF, initialized from any NMF model, where factors may be "linked" to certain samples.
#'
#' @inheritParams run_nmf
#' @param w initial matrix for 'w', usually taken from the result of \code{run_nmf}, of dimensions \code{nrow(A) x rank}.
#' @param link_h matrix giving the linkage weight (usually in the range \code{(0, 1)}) of dimensions \code{rank x ncol(A)}.
#' @param link_w matrix giving the linkage weight of dimensions \code{nrow(A) x rank}.
#'
run_linked_nmf <- function(A, w, link_h = NULL, link_w = NULL, tol = 1e-4, maxit = 100, verbose = TRUE, L1 = 0.01, L2 = 0, threads = 0) {
  if (is.null(link_h) & is.null(link_w)) {
    stop("both link_h and link_w cannot be NULL. Specify at least one linking matrix.")
  }

  if (!is.null(link_h) & nrow(link_h) != ncol(w)) {
    stop("number of rows in 'link_h' must be equal to the nubmer of columns in 'w'")
  }

  if (!is.null(link_h) & ncol(link_h) != ncol(A)) {
    stop("number of columns in 'link_h' must be equal to the number of columns in 'A'")
  }

  if (!is.null(link_w) & ncol(link_w) != ncol(w)) {
    stop("number of columns in 'link_w' must be equal to the nubmer of columns in 'w'")
  }

  if (!is.null(link_w) & nrow(link_w) != nrow(A)) {
    stop("number of rows in 'link_w' must be equal to the number of rows in 'A'")
  }

  if (is.null(link_h)) {
    link_h <- matrix(0, 1, 1)
  }

  if (is.null(link_w)) {
    link_w <- matrix(0, 1, 1)
  }

  if (L1 >= 1) {
    stop("L1 penalty must be strictly in the range (0, 1]")
  }

  if (nrow(w) != nrow(A)) {
    stop("number of rows in 'w' must be equal to the number of rows in 'A'")
  }

  link_h <- as.matrix(link_h)
  link_w <- as.matrix(link_w)
  A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  w <- t(w)

  model <- c_linked_nmf(A, t(A), tol, maxit, verbose, L1, L2, threads, w, link_h, link_w)
  sort_index <- order(model$d, decreasing = TRUE)
  model$d <- model$d[sort_index]
  model$w <- t(model$w)[, sort_index]
  model$h <- model$h[sort_index, ]
  model
}
