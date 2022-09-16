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

  A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")

  w_init <- matrix(stats::runif(nrow(A) * rank), rank, nrow(A))
  model <- c_nmf(A, t(A), tol, maxit, verbose, L1, L2, threads, w_init)
  sort_index <- order(model$d, decreasing = TRUE)
  model$d <- model$d[sort_index]
  model$w <- t(model$w)[, sort_index]
  model$h <- model$h[sort_index, ]
  rownames(model$w) <- rownames(A)
  colnames(model$h) <- colnames(A)
  colnames(model$w) <- rownames(model$h) <- paste0("NMF_", 1:ncol(model$w))
  model
}

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

  A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  At <- Matrix::t(A)

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
    model <- c_ard_nmf(A, At, tol, maxit, verbose > 1, L1, L2, threads, w_init[[rep]][1:df$k[[i]], ], abs(.Random.seed[[3 + rep]]), round(1 / test_density), tol_overfit, trace_test_mse)
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
#' @param verbose no output (0/FALSE), rank-level output (1/TRUE) and step size info (2) and individual model fitting updates (3)
#' @param learning_rate exponent on step size for automatic rank determination
#' @param tol tolerance of the final fit
#' @param cv_tol tolerance for cross-validation
#' @param k_init initial rank at which to begin search for local minimum. \code{k_init = 2} is a reasonable default, higher values can lead to swift convergence to a local minmum.
#' @import dplyr
#' @export
#'
ard_nmf <- function(A, k_init = 2, n_replicates = 3, tol = 1e-5, cv_tol = 1e-4, maxit = 100, verbose = 1, L1 = 0.01, L2 = 0, threads = 0,
                    test_density = 0.05, learning_rate = 0.8, tol_overfit = 1e-4, trace_test_mse = 5, detail_level = 1) {
  stopifnot("L1 penalty must be strictly in the range (0, 1]" = L1 < 1)
  stopifnot("'test_density' should not be greater than 0.2 or less than 0.01, as a general rule of thumb" = test_density < 0.2 & test_density > 0.01)
  A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  test_seed <- abs(.Random.seed[[3]])
  At <- Matrix::t(A)
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
      model <- c_ard_nmf(A, At, cv_tol, maxit, verbose > 2, L1, L2, threads, matrix(runif(nrow(A) * curr_rank), curr_rank, nrow(A)), test_seed + curr_rep, round(1 / test_density), tol_overfit, trace_test_mse)
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
      if(overfit_score >= tol_overfit) max_rank <- curr_rank
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
          if(higher_option >= max_rank) diff_higher <- floor(higher_option - best_rank)
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
  if (verbose > 0) cat("\nFitting final model at k =", best_rank)
  set.seed(test_seed)
  model <- c_nmf(A, At, tol, maxit, verbose > 2, L1, L2, threads, matrix(runif(nrow(A) * best_rank), best_rank, nrow(A)))
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

GetBestRank <- function(df) {
  if (ncol(df) == 5) {
    # condense to simple format by taking the last iteration in each model
    df <- as.data.frame(group_by(df, rep, k) %>% slice(which.max(iter)))
  }
  if (nrow(df) == 0) {
    return(2)
  } else if (nrow(df) == 1) {
    return(df$k[[1]])
  } else if (length(unique(df$rep)) == 1) {
    return(df$k[[which.min(df$test_error)]])
  } # get the lowest rank for each replicate, take the mean and floor it
  else {
    return(floor(mean((group_by(df, rep) %>% slice(which.min(test_error)))$k)))
  }
}

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

#' @param x the result of \code{cross_validate_nmf}
#' @param detail level of detail to plot
#' @rdname cross_validate_nmf
#' @export
#'
plot.cross_validate_nmf_data <- function(x, detail = 2, ...) {
  if (ncol(x) == 5 & detail == 1) {
    x <- as.data.frame(group_by(x, rep, k) %>% slice(which.max(iter)))
    x$iter <- NULL
  }
  if (ncol(x) < 5) {
    x$rep <- factor(x$rep)
    # simple format (detail_level = 1)
    # normalize each replicate to the same minimum
    for (rep in levels(x$rep)) {
      idx <- which(x$rep == rep)
      x$test_error[idx] <- x$test_error[idx] / min(x$test_error[idx])
    }
    best_rank <- GetBestRank(x)
    ggplot(x, aes(k, test_error, color = factor(rep))) +
      geom_point() +
      geom_line() +
      theme_classic() +
      labs(x = "factorization rank", y = "relative test set error", color = "replicate", caption = paste0("(best rank is k = ", best_rank, ")")) +
      theme(aspect.ratio = 1, plot.caption = element_text(hjust = 0.5)) +
      geom_vline(xintercept = best_rank, linetype = "dashed", color = "red") +
      scale_y_continuous(trans = "log10")
  } else {
    # detail_level = 2 format
    best_rank <- GetBestRank(x)
    if (length(unique(x$rep)) == 1) {
      ggplot(x, aes(k, test_error, color = iter, group = iter)) +
        geom_line() +
        scale_color_viridis_c(option = "B") +
        theme_classic() +
        theme(aspect.ratio = 1, plot.caption = element_text(hjust = 0.5)) +
        geom_vline(xintercept = best_rank, linetype = "dashed", color = "red") +
        scale_y_continuous(trans = "log10") +
        labs(x = "factorization rank", y = "test set error", color = "model iteration", caption = paste0("(best rank is k = ", best_rank, ")"))
    } else {
      ggplot(x, aes(k, test_error, color = iter, group = iter)) +
        geom_line() +
        scale_color_viridis_c(option = "B") +
        theme_classic() +
        theme(aspect.ratio = 1, plot.caption = element_text(hjust = 0.5)) +
        geom_vline(xintercept = best_rank, linetype = "dashed", color = "red") +
        scale_y_continuous(trans = "log10") +
        labs(x = "factorization rank", y = "test set error", color = "model iteration", caption = paste0("(best rank is k = ", best_rank, ")")) +
        facet_grid(cols = vars(rep))
    }
  }
}

#' Summarize contribution of sample groups to NMF factors
#'
#' Calculate the mean weight of samples in discrete and unique groups to each factor
#'
#' @param h matrix giving factors as rows and samples as columns
#' @param factor_data a factor of the same length as the number of columns in \code{h}
#' @param reorder sort results by proportion in each group (uses \code{hclust} if >2 groups)
#' @return \code{data.frame} of mean weights for each sample group within each factor of class \code{nmf_metadata_summary}. Use the \code{plot} method to visualize.
#' @rdname MetadataSummary
#' @export
#'
MetadataSummary <- function(h, factor_data, reorder = TRUE) {
  factor_data <- as.factor(factor_data)
  if (is.null(rownames(h))) rownames(h) <- paste0("factor", 1:nrow(h))
  m <- matrix(0, nrow(h), length(levels(factor_data)))
  rownames(m) <- rownames(h)
  colnames(m) <- levels(factor_data)
  for (j in 1:length(levels(factor_data))) {
    for (i in 1:nrow(h)) {
      m[i, j] <- mean(h[i, which(factor_data == levels(factor_data)[[j]])])
    }
  }
  m <- apply(m, 1, function(x) x / sum(x))
  if (length(levels(factor_data)) == 2) {
    m <- m[order(m[, 1], decreasing = TRUE), ]
  } else if (reorder) {
    m <- m[hclust(dist(m), method = "ward.D2")$order, hclust(dist(t(m)), method = "ward.D2")$order]
  }
  t(m)
  m <- as.data.frame(m)
  class(m) <- c("nmf_metadata_summary", "data.frame")
  m
}

#' @export
#' @rdname MetadataSummary
#' @importFrom reshape2 melt
#' @param ... not implemented
#'
plot.nmf_metadata_summary <- function(x, ...) {
  m <- reshape2::melt(as.matrix(x))
  colnames(m) <- c("group", "factor", "frac")
  ggplot(m, aes(x = factor(factor, levels = unique(factor)), y = frac, fill = group)) +
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "factor", y = "Representation in group") +
    scale_y_continuous(expand = c(0, 0))
}

#' @export
#' @rdname MetadataSummary
#' @name MetadataSummary
#'
.S3method("plot", "nmf_metadata_summary", plot.nmf_metadata_summary)

#' @export
#' @rdname MetadataSummary
#' @param x result of \code{MetadataSummary}
#' @importFrom reshape2 melt
#'
MetadataHeatmap <- function(x) {
  m <- reshape2::melt(as.matrix(x))
  colnames(m) <- c("factor", "group", "frac")
  ggplot(m, aes(x = factor(factor, levels = unique(factor)), y = group, fill = frac)) +
    geom_tile() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank()) +
    labs(x = "factor", y = "group", fill = "relative\ntotal weight") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_gradient2(low = "white", high = "red")
}

#' Load the pbmc3k dataset
#'
#' 2,700 peripheral blood mononuclear cells (PBMC) from 10x genomics taken from the "SeuratData" package
#'
#' @description
#' This dataset is adapted directly from the Satija lab "pbmc3k" dataset used in their popular tutorial on guided clustering. It is provided in this package for convenience since "SeuratData" is not available on CRAN.
#'
#' For more information, please see their documentation.
#'
#' @returns Seurat object with \code{$cell_type} info in the \code{meta.data} slot.
#'
#' @export
#'
get_pbmc3k_data <- function() {
  data(pbmc3k)
  pbmc3k
  A <- CreateSeuratObject(counts = new("dgCMatrix", i = pbmc3k$i, p = pbmc3k$p, Dim = pbmc3k$Dim, Dimnames = pbmc3k$Dimnames, x = as.numeric(inverse.rle(pbmc3k$x))))
  A@meta.data$cell_type <- pbmc3k$cell_type
  A
}

#' Compressed form of pbmc3k dataset
#'
#' @description See \code{\link{get_pbmc3k_data}}
#'
#' @md
#' @docType data
#' @usage data(pbmc3k)
#' @format compressed version of the \code{dgCMatrix}, use \code{\link{get_pbmc3k_data}} to use this dataset.
"pbmc3k"
