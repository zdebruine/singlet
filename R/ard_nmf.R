#' Automatic Rank Determination NMF
#'
#' ARD NMF quickly finds the optimal rank for an NMF model using an exponentially variable learning rate and basic coordinate descent.
#'
#' @inheritParams cross_validate_nmf
#'
#' @param verbose no output (0/FALSE), rank-level output (1/TRUE) and step size info (2) and individual model fitting updates (3)
#' @param learning_rate exponent on step size for automatic rank determination
#' @param tol tolerance of the final fit
#' @param cv_tol tolerance for cross-validation
#' @param k_init initial rank at which to begin search for local minimum. \code{k_init = 2} is a reasonable default, higher values can lead to swift convergence to a local minmum.
#' @param k_max maximum rank to consider during automatic rank determination
#' @details
#' 
#' If running ard_nmf() standalone, the following coercion can be useful:
#'
#' res <- ard_nmf(data_matrix, ...)
#' plot(res$cv_data) # rank finding
#' nmfres <- as(res, "nmf") # other
#'
#' This coercion allows AnnotateNMF, AnnotationPlot, etc. to work on `nmfres`
#' directly, rather than assuming a Seurat-like class structure is present.
#' The coercion simply checks the dimensions of res$w, res$d, and res$h,
#' then shoves all other list elements from res into nmfres@misc.
#'
#' @import dplyr
#'
#' @export
#'
ard_nmf <- function(A, k_init = 2, k_max = 100, n_replicates = 1, tol = 1e-5, cv_tol = 1e-4,
                    maxit = 100, verbose = 1, L1 = 0.01, L2 = 0, threads = 0,
                    test_density = 0.05, learning_rate = 1,
                    tol_overfit = 1e-3, trace_test_mse = 1) {
  stopifnot("L1 penalty must be strictly in the range (0, 1]" = L1 < 1)

  if(test_density > 0.02 | test_density < 0.01){
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
    w_init <- lapply(1:n_replicates, function(x) matrix(stats::runif(nrow(A[[1]]) * k_max), k_max, nrow(A[[1]])))
    sparse_list <- TRUE
    rn <- rownames(A[[1]])
    cn <- rownames(At[[1]])
  } else {
    dense_mode <- TRUE
    if (class(A)[[1]] != "matrix") {
      if (verbose > 0) cat("running with sparse optimization\n")
      A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
      At <- Matrix::t(A)
      dense_mode <- FALSE
    } else {
      if (verbose > 0) cat("running with dense optimization\n")
      At <- t(A)
    }
    w_init <- lapply(1:n_replicates, function(x) matrix(stats::runif(nrow(A) * k_max), k_max, nrow(A)))
    sparse_list <- FALSE
    rn <- rownames(A)
    cn <- colnames(A)
  }
  test_seed <- abs(.Random.seed[[3]])
  df <- data.frame("k" = integer(), "rep" = integer(), "test_error" = double(), "iter" = integer(), "tol" = double(), "overfit_score" = double())
  class(df) <- c("cross_validate_nmf_data", "data.frame")
  if (is.null(k_init) || is.na(k_init) || k_init < 2) k_init <- 2
  for (curr_rep in 1:n_replicates) {
    if (verbose >= 1 && n_replicates > 1) cat("\nREPLICATE ", curr_rep, "/", n_replicates, "\n")
    
    step_size <- 1
    curr_rank <- k_init
    while (step_size >= 1 && curr_rank < k_max) {

      # as long as we are at a rank less than the rank of the matrix and step size is 1 or greater, continue trying to find the best rank
      if (verbose > 0) cat("k =", curr_rank, ", rep =", curr_rep, "\n")
      # compute the test set reconstruction error at curr_rank
      set.seed(test_seed)
      w_init_this <- w_init[[curr_rep]][1:curr_rank, ]
      if(!sparse_list){
        if (dense_mode) {
          model <- c_ard_nmf_dense(A, At, cv_tol, maxit, verbose > 2, L1, L2, threads, w_init_this, test_seed + curr_rep, round(1 / test_density), tol_overfit, trace_test_mse)
        } else {
          model <- c_ard_nmf(A, At, cv_tol, maxit, verbose > 2, L1, L2, threads, w_init_this, test_seed + curr_rep, round(1 / test_density), tol_overfit, trace_test_mse)
        }
      } else {
        model <- c_ard_nmf_sparse_list(A, At, cv_tol, maxit, verbose > 2, L1, L2, threads, w_init_this, test_seed + curr_rep, round(1 / test_density), tol_overfit, trace_test_mse)
      }
      err <- tail(model$test_mse, n = 1L)
      overfit_score <- tail(model$score_overfit, n = 1L)
      df <- rbind(df, data.frame("k" = as.integer(curr_rank), "rep" = as.integer(curr_rep), "test_error" = model$test_mse, "iter" = model$iter, "tol" = model$tol, "overfit_score" = overfit_score))
      if (verbose > 1) {
        cat("   test_error =", sprintf(err, fmt = "%#.4e"), "\n")
        if (overfit_score == 0) cat("   not overfit\n")
        if (overfit_score > 0 & overfit_score < tol_overfit) cat("   possibly overfit (overfit_score =", sprintf(overfit_score, fmt = "%#.4e"), ")\n")
        if (overfit_score >= tol_overfit) cat("   overfit (overfit_score =", sprintf(overfit_score, fmt = "%#.4e"), ")\n")
      }
      if (overfit_score >= tol_overfit) k_max <- curr_rank
      # now decide what the next rank will be
      df_rep <- subset(df, rep == curr_rep)
      df_rep <- df_rep[order(df_rep$k), ]
      best_rank <- GetBestRank(subset(df_rep, k < k_max))
      df_rep <- group_by(df_rep, k) %>% slice(which.max(iter))
      rank_ind <- which(df_rep$k == best_rank)
      if (verbose > 1) cat("   best rank in replicate =", best_rank, "\n\n")
      if (rank_ind == nrow(df_rep)) {
        # best rank is the largest rank fit so far
        step_size <- step_size * (1 + learning_rate)
        curr_rank <- best_rank + floor(step_size)
      } else if (rank_ind == 1) {
        if (best_rank == 1) {
          cat("best rank was determined to be k = 1. You probably need more data or a better method for data normalization.\n")
          break
        }
        if (floor(step_size) < best_rank) {
          curr_rank <- best_rank - floor(step_size)
          step_size <- step_size * (learning_rate + 1)
        } else {
          curr_rank <- floor(best_rank / 2)
        }
      } else {
        next_lower_rank <- df_rep$k[[rank_ind - 1]]
        next_higher_rank <- df_rep$k[[rank_ind + 1]]
        diff_lower <- best_rank - next_lower_rank
        diff_higher <- next_higher_rank - best_rank
        higher_option <- best_rank + floor(diff_higher / 2)
        lower_option <- best_rank - floor(diff_lower / 2)
        if (diff_lower <= 1 && diff_higher <= 1) {
          break
        } else if (diff_lower >= diff_higher) {
          curr_rank <- lower_option
        } else if (diff_higher > diff_lower) {
          curr_rank <- higher_option
        }
      }
    }
  }
  df$rep <- factor(df$rep)
  class(df) <- c("cross_validate_nmf_data", "data.frame")

  best_rank <- GetBestRank(df, tol_overfit)

  # learn final nmf model
  if (verbose > 1) cat("\nUnmasking test set")
  if (verbose > 0) cat("\nFitting final model at k =", best_rank, "\n")

  set.seed(test_seed)
  w_init_this <- w_init[[1]][1:best_rank, ]
  if(!sparse_list){
    if (dense_mode) {
      model <- c_nmf_dense(A, At, tol, maxit, verbose > 2, L1, L2, threads, w_init_this)
    } else {
      model <- c_nmf(A, At, tol, maxit, verbose > 2, L1, L2, threads, w_init_this)
    }
  } else {
    model <- c_nmf_sparse_list(A, At, tol, maxit, verbose > 2, L1, L2, threads, w_init_this)
  }

  model$cv_data <- df
  sort_index <- order(model$d, decreasing = TRUE)
  model$d <- model$d[sort_index]
  model$w <- t(model$w)[, sort_index]
  model$h <- model$h[sort_index, ]
  rownames(model$w) <- rn
  colnames(model$h) <- cn
  colnames(model$w) <- rownames(model$h) <- paste0("NMF_", 1:ncol(model$w))
  model
}
