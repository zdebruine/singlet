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
#' @rdname run_nmf
#' @importFrom stats runif
#' @export
#'
run_nmf <- function(A, rank, tol = 1e-4, maxit = 100, verbose = TRUE, L1 = 0.01, L2 = 0, threads = 0) {
  if (L1 >= 1) {
    stop("L1 penalty must be strictly in the range (0, 1]")
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
    At <- lapply(1:length(A), function(i){
      do.call(rbind, lapply(1:length(A), function(j) t(A[[j]][block_sizes[i]:(block_sizes[i+1] - 1), ])))
    })
    if (verbose > 0) cat("running with sparse optimization\n")
    w_init <- matrix(stats::runif(nrow(A[[1]]) * rank), rank, nrow(A[[1]]))
    model <- c_nmf_sparse_list(A, At, tol, maxit, verbose, L1, L2, threads, w_init)
    rn <- rownames(A[[1]])
    cn <- rownames(At[[1]])
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
    
    w_init <- matrix(stats::runif(nrow(A) * rank), rank, nrow(A))
    if (dense_mode) {
      model <- c_nmf_dense(A, At, tol, maxit, verbose, L1, L2, threads, w_init)
    } else {
      model <- c_nmf(A, At, tol, maxit, verbose, L1, L2, threads, w_init)
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
