#' @rdname RunLNMF
#'
#' @export
#'
RunLNMF <- function(object, ...) {
  UseMethod("RunLNMF")
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


#' Run Linked NMF on a Seurat object
#'
#' @description Run a Linked Non-negative Matrix Factorization to separate shared and unique signals for integration or signature extraction.
#'
#' @inheritParams RunNMF.Seurat
#' @param split.by column name in \code{@meta.data} giving a \code{Factor} with multiple levels for splitting. Data will be weighted such that each group contributes equally to the LNMF model.
#' @param link.cutoff if the relative contribution of samples in any given group to a factor falls below \code{link.cutoff}, unlink it from the factor. \code{link.cutoff = 1} means a factor must contribute exactly equally before being unlinked.
#' @param link.balance after initial linking step, weight all shared factors such that each dataset is represented as equally as possible in each factor.
#' @param link.balance.tol relative change in factor representation within sample groups between balancing iterations at which to call convergence.
#' @param link.balance.rate proportion of difference between current factor weight and equal representation of factors in a sample group to target in a single iteration (default 0.1).
#' @param balance.maxit maximum number of iterations for the balancing step
#' @param reduction.use NMF reduction to use for initializing the linked factorization.
#' @param reduction.name name to store resulting DimReduc object as
#' @param reduction.key key for resulting DimReduc
#' @param verbose print fitting progress to console
#' @param tol tolerance of the fit (correlation distance of the model across consecutive iterations).
#' @param maxit maximum number of fitting iterations
#' @param L1 L1/LASSO penalty to increase sparsity of the model
#' @param L2 L2/Ridge-like penalty to increase angles between factors
#' @param threads number of threads to use (0 = let OpenMP use all available threads)
#' @param reduction reduction to use for metadata analysis
#'
#' @return a Seurat object with the NMF model stored in the reductions slot
#'
#' @details Use \code{set.seed()} to guarantee reproducibility!
#'
#' @aliases RunLNMF
#' @rdname RunLNMF
#'
#' @seealso \code{\link{RunNMF}}, \code{\link{RankPlot}}, \code{\link{MetadataSummary}}
#'
#' @export
#'
RunLNMF.Seurat <- function(object,
                           split.by,
                           reduction.use = "nmf",
                           reduction.name = "lnmf",
                           reduction.key = "LNMF_",
                           verbose = TRUE,
                           link.cutoff = 0.5,
                           tol = 1e-5,
                           maxit = 100,
                           L1 = 0.01,
                           L2 = 0,
                           threads = 0,
                           link.balance.tol = 0,
                           balance.maxit = 100,
                           link.balance.rate = 0.1, ...) {
  link_w <- NULL
  w <- object@reductions[[reduction.use]]@feature.loadings
  h <- object@reductions[[reduction.use]]@cell.embeddings
  assay <- object@reductions[[reduction.use]]@assay.used
  A <- object@assays[[assay]]@data
  rnames <- rownames(A)
  cnames <- colnames(A)
  transpose_model <- FALSE
  if (!is.null(split.by)) {
    split.by <- as.integer(as.numeric(as.factor(object@meta.data[[split.by]]))) - 1
    if (length(split.by) == nrow(A)) {
      A <- weight_by_split(t(A), split.by, length(unique(split.by)))
      transpose_model <- TRUE
    } else if (length(split.by) == ncol(A)) {
      A <- weight_by_split(A, split.by, length(unique(split.by)))
    } else {
      stop("length of 'split.by' was not equal to one of the dimensions of the input matrix")
    }
  } else {
    stop("no value specified for 'split.by'")
  }
  At <- t(A)

  # unlink factors from samples in a given group if the group representation falls below (1 - link.cutoff) / n_groups
  # matrix of groups by factors, giving mean weight of each sample in that group in that factor
  levels <- unique(split.by)
  m <- matrix(0, ncol(h), length(levels))
  for (factor in 1:nrow(m)) {
    for (level in 1:length(levels)) {
      m[factor, level] <- mean(h[which(split.by == levels[[level]]), factor])
    }
  }
  m <- t(apply(m, 1, function(x) x / sum(x))) * length(levels) < link.cutoff

  # construct linking matrix
  link_h <- matrix(1, ncol(h), nrow(h))
  for (factor in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      if (m[factor, j] == TRUE) {
        # unlink this factor from these samples
        link_h[factor, which(split.by == levels[[j]])] <- 0
      }
    }
  }

  # need to do a similar step for link_w for multi-modal integration
  link_w <- matrix(1, nrow(w), ncol(w))

  lnmf_model <- run_linked_nmf(A, w, link_h, link_w, tol, maxit, verbose, L1, L2, threads)

  # balancing step.
  # not sure if this is a good idea.
  if (link.balance.tol != 1 & link.balance.tol != 0) {
    lnmf_model$w <- t(lnmf_model$w)
    if (verbose > 0) {
      cat("balancing...\n")
    }
    # get factor representation
    m <- MetadataSummary(lnmf_model$h, split.by, FALSE)
    prev.balance.tol <- 1
    v <- as.vector(abs(1 / length(levels) - m))
    curr.balance.tol <- mean(v[v != 0 & v != 1])
    this.balance.tol <- abs(curr.balance.tol - prev.balance.tol) / (curr.balance.tol + prev.balance.tol)
    balance.iter <- 0
    while (this.balance.tol > link.balance.tol & balance.iter < balance.maxit) {
      # construct new linking matrix to correct for unequal factor representation
      for (i in 1:nrow(m)) {
        for (j in 1:ncol(m)) {
          if (m[i, j] != 1 & m[i, j] != 0) {
            if (m[i, j] < 0.5) {
              m[i, j] <- 1 + (0.5 - m[i, j]) * link.balance.rate
            } else {
              m[i, j] <- 1
            }
          }
        }
      }

      link_h <- matrix(1, ncol(h), nrow(h))
      for (factor in 1:nrow(m)) {
        for (j in 1:ncol(m)) {
          link_h[factor, which(split.by == levels[[j]])] <- m[factor, j]
        }
      }
      # run linked nmf
      lnmf_model <- c_nmf(A, At, tol, 1L, FALSE, L1, L2, threads, lnmf_model$w, link_h, 1, 0)
      m <- MetadataSummary(lnmf_model$h, split.by, FALSE)
      v <- as.vector(abs(1 / length(levels) - m))
      curr.balance.tol <- mean(v[v != 0 & v != 1])
      this.balance.tol <- abs(curr.balance.tol - prev.balance.tol) / (curr.balance.tol + prev.balance.tol)
      prev.balance.tol <- curr.balance.tol
      if (verbose > 0) cat("it: ", balance.iter, ", tol:", this.balance.tol, "\n")
      balance.iter <- balance.iter + 1
    }
    lnmf_model$w <- t(lnmf_model$w)
  }
  if (transpose_model) {
    lnmf_model <- list("w" = t(lnmf_model$h), "d" = lnmf_model$d, "h" = t(lnmf_model$h))
  }
  rownames(lnmf_model$h) <- colnames(lnmf_model$w) <- paste0(reduction.key, 1:ncol(lnmf_model$w))
  rownames(lnmf_model$w) <- rnames
  colnames(lnmf_model$h) <- cnames
  object@reductions[[reduction.name]] <- new("DimReduc",
    cell.embeddings = t(lnmf_model$h),
    feature.loadings = lnmf_model$w,
    assay.used = assay,
    stdev = as.vector(lnmf_model$d),
    key = reduction.key,
    misc = list("link_matrix" = link_h)
  )
  object
}


#' Run LNMF on a Seurat object
#'
#' S3 method for Seurat that runs the \code{singlet::RunLNMF} function.
#'
#' @method RunLNMF Seurat
#' @rdname RunLNMF
#' @name RunLNMF
#'
#' @export
#'
.S3method("RunLNMF", "Seurat", RunLNMF.Seurat)
