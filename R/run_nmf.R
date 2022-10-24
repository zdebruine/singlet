#' @title Run Non-negative Matrix Factorization
#'
#' @description Run NMF on a sparse matrix with automatic rank determination by cross-validation
#'
#' @param A sparse matrix (ideally variance-stabilized) of data for genes x cells (rows x columns)
#' @param rank factorization rank
#' @param tol tolerance of the fit (1e-5 for publication quality, 1e-3 for cross-validation)
#' @param maxit maximum number of iterations
#' @param verbose verbosity level
#' @param L1 L1/LASSO penalty to increase sparsity of model
#' @param L2 L2/Ridge penalty to increase angles between factors
#' @param threads number of threads for parallelization across CPUs, 0 = use all available threads
#' @rdname run_nmf
#' @importFrom stats runif
#' @export
#'
run_nmf <- function(A, rank, tol = 1e-4, maxit = 100, verbose = TRUE, L1 = 0.01, L2 = 0, threads = 0) {
  if (L1 >= 1) {
    stop("L1 penalty must be strictly in the range (0, 1]")
  }

  if (class(A)[[1]] != "matrix") {
    if (verbose > 0) cat("running with sparse optimization\n")
    A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    At <- Matrix::t(A)
    dense_mode <- FALSE
  } else {
    if (verbose > 0) cat("running with dense optimization\n")
    At <- t(A)
  }

  w_init <- matrix(stats::runif(nrow(A) * rank), rank, nrow(A))
  if (dense_mode) {
    model <- c_nmf_dense(A, At, tol, maxit, verbose, L1, L2, threads, w_init)
  } else {
    model <- c_nmf(A, At, tol, maxit, verbose, L1, L2, threads, w_init)
  }
  sort_index <- order(model$d, decreasing = TRUE)
  model$d <- model$d[sort_index]
  model$w <- t(model$w)[, sort_index]
  model$h <- model$h[sort_index, ]
  rownames(model$w) <- rownames(A)
  colnames(model$h) <- colnames(A)
  colnames(model$w) <- rownames(model$h) <- paste0("NMF_", 1:ncol(model$w))
  model
}
