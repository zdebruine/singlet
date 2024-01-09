#' Determine best rank for NMF using cross-validation
#'
#' @description Find the rank that minimizes the mean squared error of test set reconstruction using cross-validation.
#'
#' @inheritParams run_nmf
#' @param ranks a vector of ranks at which to fit a model and compute test set reconstruction error
#' @param n_replicates number of random test sets
#' @param test_density fraction of values to include in the test set
#' @param tol_overfit stopping criterion, maximum increase in test set reconstruction error at any iteration compared to test set reconstruction error at \code{trace_test_mse}
#' @param trace_test_mse first iteration at which to calculate test set reconstruction error, and the error to compare all later iterations to when determining whether overfitting has occurred.
#' @return a \code{data.frame} of test set reconstruction error vs. rank of class \code{nmf_cross_validate_data}. Use \code{plot} method to visualize or \code{min} to compute optimal rank.
#' @rdname cross_validate_nmf
#' @param ... additional arguments (not implemented)
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats runif
#'
cross_validate_nmf <- function(A, ranks, n_replicates = 3, tol = 1e-4, maxit = 100, verbose = 1, L1 = 0.01, L2 = 0, threads = 0, test_density = 0.05, tol_overfit = 1e-4, trace_test_mse = 5) {
  if (L1 >= 1) {
    stop("L1 penalty must be strictly in the range (0, 1]")
  }

  if (test_density > 0.2 | test_density < 0.01) {
    warning("'test_density' should not be greater than 0.2 or less than 0.01, as a general rule of thumb")
  }

  if("list" %in% class(A)){
    # check that number of rows is identical
    if(var(sapply(A, nrow)) != 0) 
      stop("number of rows in all provided 'A' matrices are not identical")
    if(!all(sapply(A, function(x) class(x) == "dgCMatrix")))
      stop("if providing a list, you must provide a list of all 'dgCMatrix' objects")
    if(!is.null(rownames(A[[1]]))){
      if(!all(sapply(A, function(x) all.equal(rownames(x), rownames(A[[1]]))))) stop("rownames of all dgCMatrix objects in list must be identical")
    }
    
    # generate a distributed transpose
    if(verbose > 0) cat("generating a distributed transpose of input matrix list\n")
    block_sizes <- floor(c(seq(1, nrow(A[[1]]), nrow(A[[1]]) /(length(A))), nrow(A[[1]]) + 1))
    At <- list()
    if(verbose > 0) pb <- txtProgressBar(min = 0, max = length(A))
    for(i in 1:length(A)){
      At[[i]] <- list()
      for(j in 1:length(A)){
        At[[i]][[j]] <- t(A[[j]][block_sizes[i]:(block_sizes[i+1] - 1), ])      
      }
      At[[i]] <- do.call(rbind, At[[i]])
      if(verbose > 0) setTxtProgressBar(pb, i)
    }
    if(verbose > 0) close(pb)
    if (verbose > 0) cat("running with sparse optimization\n")
    w_init <- lapply(1:n_replicates, function(x) matrix(stats::runif(nrow(A[[1]]) * max(ranks)), max(ranks), nrow(A[[1]])))
    sparse_list <- TRUE
  } else {
    if (class(A)[[1]] != "matrix") {
      if (verbose > 0) cat("running with sparse optimization\n")
      A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
      At <- Matrix::t(A)
      dense_mode <- FALSE
    } else {
      if (verbose > 0) cat("running with dense optimization\n")
      At <- t(A)
      dense_mode <- TRUE
    }
    w_init <- lapply(1:n_replicates, function(x) matrix(stats::runif(nrow(A) * max(ranks)), max(ranks), nrow(A)))
    sparse_list <- FALSE
  }
    
  df <- expand.grid("k" = ranks, "rep" = 1:n_replicates)
  df2 <- list()
  df$test_error <- 0
  if (verbose == 1) {
    pb <- utils::txtProgressBar(min = 0, max = nrow(df), style = 3)
  }
  for (i in 1:nrow(df)) {
    rep <- df$rep[[i]]
    if (verbose > 1) {
      cat(paste0("k = ", df$k[[i]], ", rep = ", rep, " (", i, "/", nrow(df), "):\n"))
    }
    if(!sparse_list){
      if (dense_mode) {
        model <- c_ard_nmf_dense(A, At, tol, maxit, verbose > 1, L1, L2, threads = threads, w_init[[rep]][1:df$k[[i]], ], abs(.Random.seed[[3 + rep]]), round(1 / test_density), tol_overfit, trace_test_mse)
      } else {
        model <- c_ard_nmf(A, At, tol, maxit, verbose > 1, L1, L2, threads = threads, w_init[[rep]][1:df$k[[i]], ], abs(.Random.seed[[3 + rep]]), round(1 / test_density), tol_overfit, trace_test_mse)
      }
    } else {
      model <- c_ard_nmf_sparse_list(A, At, tol, maxit, verbose > 1, L1, L2, threads = threads, w_init[[rep]][1:df$k[[i]], ], abs(.Random.seed[[3 + rep]]), round(1 / test_density), tol_overfit, trace_test_mse)
    }
    df$test_error[[i]] <- model$test_mse[[length(model$test_mse)]]
    df2[[length(df2) + 1]] <- data.frame("k" = df$k[[i]], "rep" = df$rep[[i]], "test_error" = model$test_mse, "iter" = model$iter, "tol" = model$tol)
    if (verbose == 1) utils::setTxtProgressBar(pb, i)
    if (verbose > 1) cat(paste0("test set error: ", sprintf(df$test_error[[i]], fmt = "%#.4e"), "\n\n"))
    
    if (model$test_mse[[length(model$test_mse)]] / model$test_mse[[1]] > (1 + tol_overfit)) {
      if (verbose > 1) cat(paste0("overfitting detected, lower rank recommended\n"))
    }
  }
  if (verbose == 1) close(pb)

  df$rep <- factor(df$rep)
  class(df) <- c("cross_validate_nmf_data", "data.frame")
  df2 <- do.call(rbind, df2)
  class(df2) <- c("cross_validate_nmf_data", "data.frame")
  df2
}
