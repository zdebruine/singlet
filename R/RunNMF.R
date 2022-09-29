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
  
  if(class(A)[[1]] != "matrix"){
    if(verbose > 0) cat("running with sparse optimization\n")
    A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    At <- Matrix::t(A)
    dense_mode <- FALSE
  } else {
    if(verbose > 0) cat("running with dense optimization\n")
    At <- t(A)
  }
  
  w_init <- matrix(stats::runif(nrow(A) * rank), rank, nrow(A))
  if(dense_mode){
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
#' @param detail_level \code{1} for returning best test set reconstruction error for each factorization, \code{2} to return test set reconstruction error at all \code{trace_test_mse} iterations in each factorization
#' @return a \code{data.frame} of test set reconstruction error vs. rank of class \code{nmf_cross_validate_data}. Use \code{plot} method to visualize or \code{min} to compute optimal rank.
#' @rdname cross_validate_nmf
#' @param ... additional arguments (not implemented)
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats runif
#'
cross_validate_nmf <- function(A, ranks, n_replicates = 3, tol = 1e-4, maxit = 100, verbose = 1, L1 = 0.01, L2 = 0,
                               threads = 0, test_density = 0.05, tol_overfit = 1e-4, trace_test_mse = 5, detail_level = 1) {
  if (L1 >= 1) {
    stop("L1 penalty must be strictly in the range (0, 1]")
  }

  if (test_density > 0.2 | test_density < 0.01) {
    stop("'test_density' should not be greater than 0.2 or less than 0.01, as a general rule of thumb")
  }

  if(class(A)[[1]] != "matrix"){
    if(verbose > 0) cat("running with sparse optimization\n")
    A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    At <- Matrix::t(A)
    dense_mode <- FALSE
  } else {
    if(verbose > 0) cat("running with dense optimization\n")
    At <- t(A)
  }

  df <- expand.grid("k" = ranks, "rep" = 1:n_replicates)
  df2 <- list()
  w_init <- lapply(1:n_replicates, function(x) matrix(stats::runif(nrow(A) * max(ranks)), max(ranks), nrow(A)))
  df$test_error <- 0
  if (verbose == 1) {
    pb <- utils::txtProgressBar(min = 0, max = nrow(df), style = 3)
  }
  for (i in 1:nrow(df)) {
    rep <- df$rep[[i]]
    if (verbose > 1) {
      cat(paste0("k = ", df$k[[i]], ", rep = ", rep, " (", i, "/", nrow(df), "):\n"))
    }
    if(dense_mode){
      model <- c_ard_nmf_dense(A, At, tol, maxit, verbose > 1, L1, L2, threads, w_init[[rep]][1:df$k[[i]], ], abs(.Random.seed[[3 + rep]]), round(1 / test_density), tol_overfit, trace_test_mse)
    } else {
      model <- c_ard_nmf(A, At, tol, maxit, verbose > 1, L1, L2, threads, w_init[[rep]][1:df$k[[i]], ], abs(.Random.seed[[3 + rep]]), round(1 / test_density), tol_overfit, trace_test_mse)
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
  if (detail_level == 1) {
    return(df)
  } else {
    return(df2)
  }
}

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
#'
#' @import dplyr
#'
#' @export
#'
ard_nmf <- function(A, k_init = 2, n_replicates = 3, tol = 1e-5, cv_tol = 1e-4,
                    maxit = 100, verbose = 1, L1 = 0.01, L2 = 0, threads = 0,
                    test_density = 0.05, learning_rate = 0.8,
                    tol_overfit = 1e-4, trace_test_mse = 5, detail_level = 1) {
  stopifnot("L1 penalty must be strictly in the range (0, 1]" = L1 < 1)

  stopifnot("'test_density' should not be greater than 0.2 or less than 0.01, as a general rule of thumb" = test_density < 0.2 & test_density > 0.01)
  dense_mode <- TRUE
  if(class(A)[[1]] != "matrix"){
    if(verbose > 0) cat("running with sparse optimization\n")
    A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    At <- Matrix::t(A)
    dense_mode <- FALSE
  } else {
    if(verbose > 0) cat("running with dense optimization\n")
    At <- t(A)
  }
  test_seed <- abs(.Random.seed[[3]])
  df <- data.frame("k" = integer(), "rep" = integer(), "test_error" = double())
  df2 <- list()
  if (is.null(k_init) || is.na(k_init) || k_init < 2) k_init <- 2
  for (curr_rep in 1:n_replicates) {
    if (verbose >= 1 && n_replicates > 1) cat("\nREPLICATE ", curr_rep, "/", n_replicates, "\n")
    step_size <- 1
    curr_rank <- k_init
    max_rank <- ncol(A)
    while (step_size >= 1 && curr_rank < ncol(A)) {
      # as long as we are at a rank less than the rank of the matrix and step size is 1 or greater, continue trying to find the best rank
      if (verbose > 0) cat("k =", curr_rank, ", rep =", curr_rep, "\n")
      # compute the test set reconstruction error at curr_rank
      set.seed(test_seed)
      if(dense_mode){
        model <- c_ard_nmf_dense(A, At, cv_tol, maxit, verbose > 2, L1, L2, threads, matrix(runif(nrow(A) * curr_rank), curr_rank, nrow(A)), test_seed + curr_rep, round(1 / test_density), tol_overfit, trace_test_mse)
      } else {
        model <- c_ard_nmf(A, At, cv_tol, maxit, verbose > 2, L1, L2, threads, matrix(runif(nrow(A) * curr_rank), curr_rank, nrow(A)), test_seed + curr_rep, round(1 / test_density), tol_overfit, trace_test_mse)
      }
      err <- tail(model$test_mse, n = 1L)
      overfit_score <- model$score_overfit
      df <- rbind(df, data.frame("k" = as.integer(curr_rank), "rep" = as.integer(curr_rep), "test_error" = err))
      df2[[length(df2) + 1]] <- data.frame("k" = as.integer(curr_rank), "rep" = as.integer(curr_rep), "test_error" = model$test_mse, "iter" = model$iter, "tol" = model$tol)
      if (verbose > 1) {
        cat("   test_error =", sprintf(err, fmt = "%#.4e"), "\n")
        if (overfit_score == 0) cat("   not overfit\n")
        if (overfit_score > 0 & overfit_score < tol_overfit) cat("   possibly overfit (overfit_score =", sprintf(overfit_score, fmt = "%#.4e"), ")\n")
        if (overfit_score >= tol_overfit) cat("   overfit (overfit_score =", sprintf(overfit_score, fmt = "%#.4e"), ")\n")
      }
      if (overfit_score >= tol_overfit) max_rank <- curr_rank
      # now decide what the next rank will be
      df_rep <- subset(df, rep == curr_rep)
      df_rep <- df_rep[order(df_rep$k), ]
      best_rank <- GetBestRank(subset(df_rep, k < max_rank))
      best_rank_pos <- which(df_rep$k == best_rank)
      if (verbose > 1) cat("   best rank in replicate =", best_rank, "\n\n")
      if (best_rank_pos == length(df_rep$k)) {
        # best rank is the largest rank fit so far
        step_size <- step_size * (1 + learning_rate)
        curr_rank <- best_rank + floor(step_size)
      } else {
        if (best_rank_pos == 1) {
          if (best_rank == 1) break
          if (floor(step_size) < best_rank) {
            curr_rank <- best_rank - floor(step_size)
            step_size <- step_size * (learning_rate + 1)
          } else {
            curr_rank <- floor(best_rank / 2)
          }
        } else {
          next_lower_rank <- df_rep$k[[best_rank_pos - 1]]
          next_higher_rank <- df_rep$k[[best_rank_pos + 1]]
          diff_lower <- best_rank - next_lower_rank
          diff_higher <- next_higher_rank - best_rank
          higher_option <- best_rank + floor(diff_higher / 2)
          if (higher_option >= max_rank) diff_higher <- floor(higher_option - best_rank)
          if (diff_lower <= 1 && diff_higher <= 1) {
            break
          } else if (diff_lower > diff_higher) {
            curr_rank <- best_rank - floor(diff_lower / 2)
          } else if (diff_higher >= diff_lower) {
            curr_rank <- best_rank + floor(diff_higher / 2)
          }
        }
      }
    }
  }
  df <- arrange(df, rep, k)
  df$rep <- factor(df$rep)
  class(df) <- c("cross_validate_nmf_data", "data.frame")

  best_rank <- GetBestRank(df)

  df2 <- do.call(rbind, df2)
  class(df2) <- c("cross_validate_nmf_data", "data.frame")

  # learn final nmf model
  if (verbose > 1) cat("\nUnmasking test set")
  if (verbose > 0) cat("\nFitting final model at k =", best_rank, "\n")
  set.seed(test_seed)
  if(dense_mode){
    model <- c_nmf_dense(A, At, tol, maxit, verbose > 2, L1, L2, threads, matrix(runif(nrow(A) * best_rank), best_rank, nrow(A)))
  } else {
    model <- c_nmf(A, At, tol, maxit, verbose > 2, L1, L2, threads, matrix(runif(nrow(A) * best_rank), best_rank, nrow(A)))
  }
  ifelse(detail_level == 1, model$cv_data <- df, model$cv_data <- df2)
  sort_index <- order(model$d, decreasing = TRUE)
  model$d <- model$d[sort_index]
  model$w <- t(model$w)[, sort_index]
  model$h <- model$h[sort_index, ]
  rownames(model$w) <- rownames(A)
  colnames(model$h) <- colnames(A)
  colnames(model$w) <- rownames(model$h) <- paste0("NMF_", 1:ncol(model$w))
  model
}



#' Run NMF on a Seurat object
#'
#' @description Run Non-negative Matrix Factorization with rank determined by CV
#'
#' @param object A Seurat object
#' @param assay Assay to use, defaults to the default assay of the first object
#' @param reduction.name Name to store resulting DimReduc object as
#' @param reduction.key Key for resulting DimReduc
#' @param verbose Level of console output (0/FALSE, 1/TRUE, 2)
#' @param k an initial rank at which to start automatic cross-validation, or a vector of ranks at which to fit.
#' @param reps number of replicates for cross-validation
#' @param test.set.density approximate density of the test set (default 0.05)
#' @param tol tolerance of the fit (correlation distance of the model across consecutive iterations). Cross-validation fits are 10x coarser than this tolerance.
#' @param maxit maximum number of fitting iterations
#' @param L1 L1/LASSO penalty to increase sparsity of the model
#' @param L2 L2/Ridge-like penalty to increase angles between factors
#' @param threads number of threads to use (0 = let OpenMP use all available threads)
#' @param split.by column name in \code{@meta.data} giving a \code{Factor} with multiple levels for splitting. Data will be weighted such that each group contributes equally to the LNMF model.
#' @param learning.rate exponent on step size for automatic rank determination
#' @param tol.overfit tolerance for increase in test set reconstruction error relative to minimum observed value during factorization
#' @param trace.test.mse during automatic rank determination, calculate test set reconstruction error every trace iterations
#' @param ... not implemented
#'
#' @return Returns a Seurat object with the NMF model stored in the reductions slot
#'
#' @details Use \code{set.seed()} to guarantee reproducibility!
#' @rdname RunNMF
#' @aliases RunNMF.Seurat
#' @name RunNMF.Seurat
#'
#' @examples
#' get_pbmc3k_data() %>%
#'   NormalizeData() %>%
#'   RunNMF() -> pbmc3k
#'
#' @seealso \code{\link{RunLNMF}}, \code{\link{RankPlot}}, \code{\link{MetadataSummary}}
#'
#' @export
#'
RunNMF.Seurat <- function(object,
                          split.by = NULL,
                          k = 2,
                          assay = NULL,
                          reps = 3,
                          tol = 1e-5,
                          L1 = 0.01,
                          L2 = 0,
                          verbose = 2,
                          reduction.name = "nmf",
                          reduction.key = "NMF_",
                          maxit = 100,
                          test.set.density = 0.05,
                          learning.rate = 0.8,
                          tol.overfit = 1e-4,
                          trace.test.mse = 5,
                          threads = 0,
                          ...) {
  if (length(k) <= 1) {
    stopifnot("k must be an integer or vector of integers" = !is.na(k) || !is.null(k))
  }

  if (is.null(assay)) {
    assay <- names(object@assays)[[1]]
  }

  # check if data has been normalized
  v <- object@assays[[assay]]@data@x
  if (sum(as.integer(v)) == sum(v)) {
    object <- PreprocessData(object, assay = assay)
  }
  A <- object@assays[[assay]]@data
  rnames <- rownames(A)
  cnames <- colnames(A)

  if (!is.null(split.by)) {
    split.by <- as.integer(as.numeric(as.factor(object@meta.data[[split.by]]))) - 1
    if (any(sapply(split.by, is.na))) {
      stop("'split.by' cannot contain NA values")
    }
    A <- weight_by_split(A, split.by, length(unique(split.by)))
  }
  At <- Matrix::t(A)
  seed.use <- abs(.Random.seed[[3]])

  if (length(k) > 1) {
    # run cross-validation at specified ranks
    cv_data <- data.frame()
    cv_data <- cross_validate_nmf(A, k, reps, tol * 10, maxit, verbose, L1, L2, threads, test.set.density, tol.overfit, trace.test.mse, 2)
    best_rank <- GetBestRank(cv_data)
    if (verbose >= 1) {
      cat("best rank: ", best_rank, "\n")
    }
    cat("\nfitting final model:\n")
    nmf_model <- run_nmf(A, best_rank, tol, maxit, verbose > 1, L1, L2, threads)
  } else if (length(k) == 1) {
    # run automatic rank determination cross-validation
    nmf_model <- ard_nmf(
      A = A,
      k_init = k,
      n_replicates = reps,
      tol = tol,
      maxit = maxit,
      verbos = verbose,
      L1 = L1,
      L2 = L2,
      threads = threads,
      test_density = test.set.density,
      learning_rate = learning.rate,
      tol_overfit = tol.overfit,
      trace_test_mse = trace.test.mse,
      detail_level = 2
    )
    cv_data <- nmf_model$cv_data
  } else {
    stop("value for 'k' was invalid")
  }
  rownames(nmf_model$h) <- colnames(nmf_model$w) <- paste0(reduction.key, 1:nrow(nmf_model$h))
  rownames(nmf_model$w) <- rnames
  colnames(nmf_model$h) <- cnames
  object@reductions[[reduction.name]] <- new("DimReduc",
    cell.embeddings = t(nmf_model$h),
    feature.loadings = nmf_model$w,
    assay.used = assay,
    stdev = nmf_model$d,
    key = reduction.key,
    misc = list("cv_data" = cv_data)
  )

  object
}


#' @rdname RunNMF
#'
#' @name RunNMF
#'
#' @export
#'
RunNMF <- function(object, ...) {
  UseMethod("RunNMF")
}


#' @rdname RunNMF
#'
#' @name RunNMF
#'
#' @export
#'
.S3method("RunNMF", "Seurat", RunNMF.Seurat)
