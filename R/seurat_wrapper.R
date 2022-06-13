#' Run NMF on a Seurat object
#' 
#' @description Run a Non-negative Matrix Factorization dimensionality reduction with rank determination by cross-validation.
#' 
#' @param object A Seurat object
#' @param assay Assay to use, defaults to the default assay of the first object
#' @param reduction.name Name to store resulting DimReduc object as
#' @param reduction.key Key for resulting DimReduc
#' @param verbose Level of console output (0/FALSE, 1/TRUE, 2)
#' @param k vector of ranks at which to fit, witholding a test set
#' @param n_replicates number of replicates for cross-validation
#' @param test_set approximate density of the test set
#' @param tol tolerance of the fit (correlation distance of the model across consecutive iterations)
#' @param maxit maximum number of fitting iterations
#' @param L1 L1/LASSO penalty to increase sparsity of the model
#' @param threads number of threads to use (0 = let OpenMP use all available threads)
#' @param split.by Attribute in \code{@metadata} for splitting, if applicable. Data will be weighted such that each group contributes equally to the NMF model.
#' 
#' @details Use \code{set.seed()} to guarantee reproducibility!
#' 
#' @return Returns a Seurat object with the NMF model stored in the reductions slot
#' 
RunNMF <- function(
  object, 
  assay = NULL, 
  reduction.name = "nmf", 
  reduction.key = "NMF_", 
  verbose = TRUE, 
  k = seq(2, 40, 2),
  n_replicates = 1,
  test.set.density = 0.05,
  tol = 1e-4,
  maxit = 100,
  L1 = 0.01,
  threads = 0,
  split.by = NULL,
  precision = "float"){
  
  if(is.null(assay))
    assay <- names(object@assays)[[1]]
  
  # check if data has been normalized
  v <- ifnb@assays[[assay]]@data@x
  if(sum(as.integer(v)) == sum(v)){
    object <- NormalizeData(object, assay = assay)
  }
  A <- object@assays[[assay]]@data
  rnames <- rownames(A)
  cnames <- colnames(A)
  
  if(!is.null(split.by)){
    split.by <- as.integer(as.numeric(as.factor(object@meta.data[[split.by]]))) - 1
    A <- weight_by_split(A, split.by, length(unique(split.by)))
  }
  At <- t(A)
  seed.use <- abs(.Random.seed[[3]])

  cv_data <- data.frame()
  if(length(k) > 1){
    cv_data <- cross_validate_nmf(A, k, n_replicates, tol, maxit, verbose, L1, threads, test.set.density, precision)
    best_rank <- GetBestRank(cv_data)
    if(verbose >= 1)
      cat("best rank: ", best_rank, "\n")
    cat("\nfitting final model:\n")
    k <- best_rank
  }
  nmf_model <- run_nmf(A, k, tol, maxit, verbose > 1, L1, threads, precision)
  rownames(nmf_model$h) <- colnames(nmf_model$w) <- paste0(reduction.key, 1:k)
  rownames(nmf_model$w) <- rnames
  colnames(nmf_model$h) <- cnames
  object@reductions[[reduction.name]] <- new("DimReduc", 
                                             cell.embeddings = t(nmf_model$h), 
                                             feature.loadings = nmf_model$w, 
                                             assay.used = assay, 
                                             stdev = nmf_model$d, 
                                             key = reduction.key, 
                                             misc = list("cv_data" = cv_data))
  
  # get model with minimum reconstruction error
  object
}

#' Run Linked NMF on a Seurat object
#' 
#' @description Run a Linked Non-negative Matrix Factorization to separate shared and unique signals for integration or signature extraction.
#' 
#' @param object Seurat object
#' @param split.by attribute in \code{@metadata} for splitting. Data will be weighted such that each group contributes equally to the LNMF model.
#' @param link.cutoff if the relative contribution of samples in any given group to a factor falls below \code{link.cutoff}, unlink it from the factor. \code{link.cutoff = 1} means a factor must contribute exactly equally before being unlinked.
#' @param link.balance weight all shared factors such that each dataset is represented equally in each factor.
#' @param reduction.use NMF reduction to use for initializing the linked factorization.
#' @param reduction.name name to store resulting DimReduc object as
#' @param reduction.key key for resulting DimReduc
#' @param verbose print fitting progress to console
#' @param tol tolerance of the fit (correlation distance of the model across consecutive iterations)
#' @param maxit maximum number of fitting iterations
#' @param L1 L1/LASSO penalty to increase sparsity of the model
#' @param threads number of threads to use (0 = let OpenMP use all available threads)
#' 
#' @details Use \code{set.seed()} to guarantee reproducibility!
#' 
#' @return Returns a Seurat object with the NMF model stored in the reductions slot
#' 
RunLNMF <- function(
  object,
  split.by,
  reduction.use = "nmf",
  reduction.name = "lnmf",
  reduction.key = "LNMF_",
  verbose = TRUE,
  link.cutoff = 0.5,
  tol = 1e-4,
  maxit = 100,
  L1 = 0.01,
  threads = 0,
  precision = "float",
  link.balance.tol = 1e-5,
  balance.maxit = 100,
  link.balance.rate = 0.1){
  
  link_w <- NULL
  w <- object@reductions[[reduction.use]]@feature.loadings
  h <- object@reductions[[reduction.use]]@cell.embeddings
  assay <- object@reductions[[reduction.use]]@assay.used
  A <- object@assays[[assay]]@data
  rnames <- rownames(A)
  cnames <- colnames(A)
  transpose_model <- FALSE
  if(!is.null(split.by)){
    split.by <- as.integer(as.numeric(as.factor(object@meta.data[[split.by]]))) - 1
    if(length(split.by) == nrow(A)){
      A <- weight_by_split(t(A), split.by. length(unique(split.by)))
      transpose_model <- TRUE
    } else if(length(split.by) == ncol(A)){
      A <- weight_by_split(A, split.by, length(unique(split.by)))
    } else stop("length of 'split.by' was not equal to one of the dimensions of the input matrix")
  } else stop("no value specified for 'split.by'")
  At <- t(A)

  # unlink factors from samples in a given group if the group representation falls below (1 - link.cutoff) / n_groups
  levels <- unique(split.by)
  m <- matrix(0, length(split.by), length(levels))
  for(j in 1:nrow(m)){
    for(i in 1:length(levels)){
      m[j, i] <- mean(h[which(split.by == levels[[i]]), j])
    }
  }
  m <- t(apply(m, 1, function(x) x / sum(x))) * length(levels) < link.cutoff

  link_h <- matrix(1, ncol(h), nrow(h))
  for(factor in 1:nrow(m)){
    for(j in 1:ncol(m)){
      if(m[factor, j] == TRUE){
        # unlink this factor from these samples
        link_h[factor, which(split.by == levels[[j]])] <- 0
      }
    }
  }

  lnmf_model <- run_linked_nmf(A, w, link_h, link_w, tol, maxit, verbose, L1, threads, precision)
  if(link.balance.tol != 1 & link.balance.tol != 0 & !is.na(link.balance.tol)){
    lnmf_model$w <- t(lnmf_model$w)
    if(verbose > 0)
      cat("balancing...\n")
    # get factor representation
    m <- MetadataSummary(lnmf_model$h, split.by, cluster = FALSE)
    prev.balance.tol <- 1
    v <- as.vector(abs(1 / length(levels) - m))
    curr.balance.tol <- mean(v[v != 0 & v != 1])
    this.balance.tol <- abs(curr.balance.tol - prev.balance.tol) / (curr.balance.tol + prev.balance.tol)
    balance.iter <- 0
    while(this.balance.tol > link.balance.tol & balance.iter < balance.maxit){
      # construct new linking matrix to correct for unequal factor representation
      for(i in 1:nrow(m)){
        for(j in 1:ncol(m)){
          if(m[i, j] != 1 & m[i, j] != 0){
            if(m[i, j] < 0.5){
              m[i, j] <- 1 + (0.5 - m[i, j]) * link.balance.rate
            } else {
              m[i, j] <- 1
            }
          }
        }
      }

      link_h <- matrix(1, ncol(h), nrow(h))
      for(factor in 1:nrow(m)){
        for(j in 1:ncol(m)){
          link_h[factor, which(split.by == levels[[j]])] <- m[factor, j]
        }
      }
      # run linked nmf
      lnmf_model <- c_nmf(A, At, tol, 1L, FALSE, L1, threads, lnmf_model$w, link_h, 1, 0, precision == "float")
      m <- MetadataSummary(lnmf_model$h, split.by, cluster = FALSE)
      v <- as.vector(abs(1 / length(levels) - m))
      curr.balance.tol <- mean(v[v != 0 & v != 1])
      this.balance.tol <- abs(curr.balance.tol - prev.balance.tol) / (curr.balance.tol + prev.balance.tol)
      prev.balance.tol <- curr.balance.tol
      if(verbose > 0) cat("it: ", balance.iter, ", tol:", this.balance.tol, "\n")
      balance.iter <- balance.iter + 1;
    }
    lnmf_model$w <- t(lnmf_model$w)
  }
  if(transpose_model){
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
                                             misc = list("link_matrix" = link_h))
  object
}

RankPlot.Seurat <- function(x, reduction = "nmf"){
  RankPlot(x@reductions[[reduction]]@misc$cv_data)
}

MetadataPlot.Seurat <- function(x, split.by, reduction = "nmf"){
  MetadataPlot(t(x@reductions[[reduction]]@cell.embeddings), x@meta.data[[split.by]])
}

GetSharedFactors <- function(x, split.by, reduction = "lnmf"){
  which(!(colnames(x@reductions[[reduction]]@cell.embeddings) %in% names(which(apply(MetadataSummary(t(x@reductions[[reduction]]@cell.embeddings), x@meta.data[[split.by]]), 1, function(x) min(x) == 0)))))
}

GetUniqueFactors <- function(x, split.by, reduction = "lnmf"){
  which((colnames(x@reductions[[reduction]]@cell.embeddings) %in% names(which(apply(MetadataSummary(t(x@reductions[[reduction]]@cell.embeddings), x@meta.data[[split.by]]), 1, function(x) min(x) == 0)))))
}
