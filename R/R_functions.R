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
  model <- c_nmf(A, t(A), tol, maxit, verbose, L1, L2, threads, w_init, matrix(), matrix(), 1, 0)
  sort_index <- order(model$d, decreasing = TRUE)
  model$d <- model$d[sort_index]
  model$w <- t(model$w)[, sort_index]
  model$h <- model$h[sort_index, ]
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
#' @return a \code{data.frame} of test set reconstruction error vs. rank of class \code{nmf_cross_validate_data}. Use \code{plot} method to visualize or \code{min} to compute optimal rank.
#' @rdname cross_validate_nmf
#' @param ... additional arguments (not implemented)
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats runif
#'
cross_validate_nmf <- function(A, ranks, n_replicates = 3, tol = 1e-4, maxit = 100, verbose = 1, L1 = 0.01, L2 = 0, threads = 0, test_density = 0.05) {
  if (L1 >= 1) {
    stop("L1 penalty must be strictly in the range (0, 1]")
  }

  if (test_density > 0.2 | test_density < 0.01) {
    stop("'test_density' should not be greater than 0.2 or less than 0.01, as a general rule of thumb")
  }

  A <- as(as(as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  At <- Matrix::t(A)

  df <- expand.grid("k" = ranks, "rep" = 1:n_replicates)
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
    df$test_error[[i]] <- c_nmf(A, At, tol, maxit, verbose > 1, L1, L2, threads, w_init[[rep]][1:df$k[[i]], ], matrix(0, 1, 1), matrix(0, 1, 1), abs(.Random.seed[[3 + rep]]), round(1 / test_density))$test_mse
    if (verbose == 1) utils::setTxtProgressBar(pb, i)
    if (verbose > 1) cat(paste0("test set error: ", sprintf(df$test_error[[i]], fmt = "%#.4e"), "\n\n"))
  }
  if (verbose == 1) close(pb)
  df$rep <- factor(df$rep)
  class(df) <- c("cross_validate_nmf_data", "data.frame")
  df
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

  model <- c_nmf(A, t(A), tol, maxit, verbose, L1, L2, threads, w, link_h, link_w, 1, 0)
  sort_index <- order(model$d, decreasing = TRUE)
  model$d <- model$d[sort_index]
  model$w <- t(model$w)[, sort_index]
  model$h <- model$h[sort_index, ]
  model
}

GetBestRank <- function(df, ...) {
  if (!("cross_validate_nmf_data" %in% class(df))) {
    stop("input data.frame must be a result of 'cross_validate_nmf'")
  }

  (df %>% group_by(k) %>% dplyr::summarize(Mean = mean(test_error)) %>% slice(which.min(Mean)))$k
}

#' @param x the result of \code{cross_validate_nmf}
#' @rdname cross_validate_nmf
#' @export
#'
plot.cross_validate_nmf_data <- function(x, ...) {

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
    geom_vline(xintercept = best_rank, linetype = "dashed", color = "red")
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
