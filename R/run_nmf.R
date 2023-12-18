#' @title Run Non-negative Matrix Factorization
#'
#' @description Run NMF on a sparse matrix with automatic rank determination by cross-validation
#'
#' @param A sparse matrix giving normalized counts for genes x cells (rows x columns), or a list of sparse matrices with equal number of rows and identical rownames
#' @param rank factorization rank
#' @param tol tolerance of the fit (1e-5 for publication quality, 1e-4 for cross-validation)
#' @param maxit maximum number of iterations
#' @param verbose verbosity level
#' @param L1 L1/LASSO penalty to increase sparsity of model
#' @param L2 L2/Ridge penalty to increase angles between factors
#' @param threads number of threads for parallelization across CPUs, 0 = use all available threads
#' @param compression_level either 2 or 3, for VCSC or IVCSC, respectively. For development purposes.
#' @rdname run_nmf
#' @importFrom stats runif
#' @export
#'
run_nmf <- function(A, rank, tol = 1e-4, maxit = 100, verbose = TRUE, L1 = 0.01, L2 = 0, threads = 0, compression_level = 3) {
  use_vcsc <- compression_level == 2

  if ("list" %in% class(A)) {
    # check that number of rows is identical
    if (var(sapply(A, nrow)) != 0) {
      stop("number of rows in all provided 'A' matrices are not identical")
    }
    if (!all(sapply(A, function(x) class(x) == "dgCMatrix"))) {
      stop("if providing a list, you must provide a list of all 'dgCMatrix' objects")
    }
    if (!is.null(rownames(A[[1]]))) {
      if (!all(sapply(A, function(x) all.equal(rownames(x), rownames(A[[1]]))))) stop("rownames of all dgCMatrix objects in list must be identical")
    }
    w_init <- matrix(stats::runif(nrow(A[[1]]) * rank), rank, nrow(A[[1]]))
    model <- run_nmf_on_sparsematrix_list(A, tol, maxit, verbose, threads, w_init, use_vcsc)
    rn <- rownames(A[[1]])
    cn <- do.call(c, lapply(A, colnames))
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
    
    if(length(L1) != 2){
      L1 <- c(L1[[1]], L1[[1]])
    }
    if(length(L2) != 2){
      L2 <- c(L2[[1]], L2[[1]])
    }

    w_init <- matrix(stats::runif(nrow(A) * rank), rank, nrow(A))
    if (dense_mode) {
      model <- c_nmf_dense(A, At, tol, maxit, verbose, L1[[1]], L1[[2]], L2[[1]], L2[[2]], threads, w_init)
    } else {
      model <- c_nmf(A, At, tol, maxit, verbose, L1[[1]], L1[[2]], L2[[1]], L2[[2]], threads, w_init)
    }
    rn <- rownames(A)
    cn <- colnames(A)
  }

  sort_index <- order(model$d, decreasing = TRUE)
  model$d <- model$d[sort_index]
  model$w <- t(model$w)[, sort_index]
  model$h <- model$h[sort_index, ]
  rownames(model$w) <- rn
  colnames(model$h) <- cn
  colnames(model$w) <- rownames(model$h) <- paste0("NMF_", 1:ncol(model$w))
  model
}

distributed_transpose <- function(A){
  library(Matrix)
  setwd("/active/debruinz_project/debruinz/CellCensusNMF")
  A <- lapply(paste0("../../CellCensus/R/chunk", 1:100, "_counts.rds"), readRDS)
  block_sizes <- floor(c(seq(1, nrow(A[[1]]), nrow(A[[1]]) / (length(A))), nrow(A[[1]]) + 1))
  for (i in 1:length(block_sizes)) {
    cat("CHUNK", i, "/100\n")
    At <- list()
    pb <- txtProgressBar(min = 0, max = length(A), style = 3)
    for (j in 1:length(A)) {
      At[[j]] <- t(A[[j]][block_sizes[i]:(block_sizes[i + 1] - 1), ])
      setTxtProgressBar(pb, j)
    }
    cat("   rbinding\n")
    At <- do.call(rbind, At)
    cat("   saving\n")
    saveRDS(At, paste0("chunk", i, "_transpose_counts.rds"))
  }
}

split_into_chunks <- function(A, n_chunks){
  breakpoints <- seq(1, ncol(A), floor(ncol(A) / n_chunks))
  breakpoints[length(breakpoints) + 1] <- ncol(A)
  result <- list()
  for(i in 1:n_chunks){
    result[[i]] <- A[,breakpoints[i]:breakpoints[i + 1]]
  }
  result
}