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
#' get_pbmc3k_data() %>% NormalizeData() %>% RunNMF() -> pbmc3k
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
    nmf_model <- ard_nmf(A=A, 
                         k_init=k, 
                         n_replicates=reps, 
                         tol=tol, 
                         maxit=maxit, 
                         verbos=verbose, 
                         L1=L1, 
                         L2=L2, 
                         threads=threads,
                         test_density=test.set.density,
                         learning_rate=learning.rate, 
                         tol_overfit=tol.overfit, 
                         trace_test_mse=trace.test.mse,
                         detail_level=2)
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
