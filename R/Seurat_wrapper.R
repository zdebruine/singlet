#' Run NMF on a Seurat object
#'
#' @description Run a Non-negative Matrix Factorization dimensionality reduction with rank determination by cross-validation.
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
#' @details Use \code{set.seed()} to guarantee reproducibility!
#' @export
#' @aliases RunNMF
#' @seealso \code{\link{RunLNMF}}, \code{\link{RankPlot}}, \code{\link{MetadataSummary}}
#' @rdname RunNMF
#' @return Returns a Seurat object with the NMF model stored in the reductions slot
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
    nmf_model <- ard_nmf(A, k, reps, tol, maxit, verbose, L1, L2, threads, test.set.density, learning.rate, tol.overfit, trace.test.mse, 2)
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

#' Normalize count data
#' 
#' Standard log-normalization equivalent to \code{Seurat::LogNormalize}, but faster and more memory-efficient
#'  
#' @param object Seurat object
#' @param assay assay in which the counts matrix resides
#' @param scale.factor value by which to multiply all columns after unit normalization and before \code{log1p} transformation
#' @param ... not implemented
#' @export
#' @rdname PreprocessData
#'
PreprocessData.Seurat <- function(object, scale.factor = 10000, assay = NULL, ...){
  if(is.null(assay)) assay <- names(object@assays)[[1]]
  if(is.null(object@assays[[assay]]@key))
    object@assays[[assay]]@key <- paste0(assay, "_")
  object@assays[[assay]] <- PreprocessData(object@assays[[assay]], ...)
  object
}

#' @rdname PreprocessData
#' @export
PreprocessData.Assay <- function(object, scale.factor = 10000, ...){
  if(ncol(object@counts) == 0){
    object@data <- PreprocessData(object@data, ...)
  } else {
    object@data <- PreprocessData(object@counts, ...)
  }
  object
}

#' @rdname PreprocessData
#' @export
PreprocessData.dgCMatrix <- function(object, scale.factor = 10000, ...){
  log_normalize(object, scale.factor, 0)
}

#' @export
#' @rdname PreprocessData
#'
PreprocessData <- function(object, scale.factor, ...) {
  UseMethod("PreprocessData")
}

#' @export
#' @rdname PreprocessData
#' @name PreprocessData
#'
.S3method("PreprocessData", "dgCMatrix", PreprocessData.dgCMatrix)

#' @export
#' @rdname PreprocessData
#' @name PreprocessData
#'
.S3method("PreprocessData", "Assay", PreprocessData.Assay)


#' @export
#' @rdname PreprocessData
#' @name PreprocessData
#'
.S3method("PreprocessData", "Seurat", PreprocessData.Seurat)


#' Project data onto a factor model
#'
#' @description Non-negative Least Squares (NNLS) projection of assay data onto a factor model for transfer learning
#'
#' @inheritParams RunNMF.Seurat
#' @param w a factor model of the same number of rows as the assay to be projected in the current object, where columns correspond to factors
#' @param ... not implemented
#' @details Use \code{set.seed()} to guarantee reproducibility!
#' @export
#' @aliases ProjectData
#' @seealso \code{\link{RunLNMF}}, \code{\link{MetadataSummary}}
#' @rdname ProjectData
#' @return Returns a Seurat object with the projection stored in the reductions slot
#'
ProjectData.Seurat <- function(object, w,
                          split.by = NULL,
                          assay = NULL,
                          L1 = 0.01,
                          L2 = 0,
                          reduction.name = "nmf_projection",
                          reduction.key = "NNLS_",
                          threads = 0,
                          ...) {
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
  
  w <- w[which(rownames(w) %in% rownames(A)), ]
  A <- A[which(rownames(A) %in% rownames(w)), ]
  
  nmf_model <- project_model(A, w, L1, L2, threads)
  nmf_model$w <- w
  rownames(nmf_model$h) <- colnames(nmf_model$w) <- paste0(reduction.key, 1:nrow(nmf_model$h))
  rownames(nmf_model$w) <- rownames(A)
  colnames(nmf_model$h) <- cnames
  idx <- order(nmf_model$d, decreasing = TRUE)
  nmf_model$h <- nmf_model$h[idx, ]
  nmf_model$d <- nmf_model$d[idx]
  nmf_model$w <- nmf_model$w[, idx]
  object@reductions[[reduction.name]] <- new("DimReduc",
                                             cell.embeddings = t(nmf_model$h),
                                             feature.loadings.projected = nmf_model$w,
                                             assay.used = assay,
                                             stdev = nmf_model$d,
                                             key = reduction.key)
  object
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
#' @seealso \code{\link{RunNMF}}, \code{\link{RankPlot}}, \code{\link{MetadataSummary}}
#' @export
#' @details Use \code{set.seed()} to guarantee reproducibility!
#' @aliases RunLNMF
#' @rdname RunLNMF
#' @return Returns a Seurat object with the NMF model stored in the reductions slot
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

#' Annotate NMF model with cell metadata
#' 
#' @details Maps factor information from the \code{meta.data} slot of a Seurat object or a specified factor to each NMF factor, and computes summary statistics
#' 
#' @inheritParams RunNMF
#' @param fields if a Seurat object is specified, fields in the \code{meta.data} slot may be specified for which to compute summary statistics. If \code{fields = NULL}, all fields that are a \code{factor} with more than one level are considered.
#' @param meta.data if a DimReduc object is specified, provide a named list containing factor(s) for which to calculate summary statistics against the NMF model
#' @export
#' @rdname AnnotateNMF
#' @aliases AnnotateNMF
#'
AnnotateNMF.DimReduc <- function(object, meta.data, ...){

  # get all factors with multiple levels from meta.data data.frame
  fields <- c()
  for(i in 1:length(names(meta.data))){
    if(is.factor(meta.data[, i]) & length(levels(meta.data[, i])) > 1){
      fields <- c(fields, names(meta.data)[[i]])
    }
  }
  if(length(fields) == 0)
    stop("there were no factors in the meta.data slot with more than 1 level")
  
  # loop through each of these factors and run tests on a grid of all levels vs. all factors
  annotations <- list()
  for(field in fields){
    v <- meta.data[[field]]
    df <- expand.grid("group" = levels(v), "factor" = 1:ncol(object@cell.embeddings))
    df$p <- 1
    df$fc <- 0
    for(i in 1:nrow(df)){
      vals <- object@cell.embeddings[, df$factor[[i]]]
      in_group <- vals[which(v == df$group[[i]])]
      out_group <- vals[-which(v == df$group[[i]])]
      fc <- mean(in_group) / mean(out_group)
      df$fc[[i]] <- ifelse(fc > 1, fc - 1, -1 / fc + 1)
      if(df$fc[[i]] > 0)
        df$p[[i]] <- t.test(in_group, out_group, "greater")$p.value
    }
    annotations[[field]] <- df
  }
  object@misc$annotations <- annotations
  object
}

#' @export
#' @rdname AnnotateNMF
#' @param reduction the reductions slot in the Seurat object containing the model to annotate
#' @aliases AnnotateNMF
#' 
AnnotateNMF.Seurat <- function(object, fields = NULL, reduction = "nmf", ...){
  if(!is.null(fields)){
    for(i in 1:length(fields)){
      if(!is.factor(object@meta.data[[fields[[i]]]])){
        stop(fields[[i]], "is not a factor!")
      } else if(length(levels(object@meta.data[[fields[[i]]]]))){
        stop(fields[[i]], "is a factor, but only has one level")
      }
    }
    object@reductions[[reduction]] <- AnnotateNMF.DimReduc(object@reductions[[reduction]], object@meta.data[, fields], ...)
  } else {
    object@reductions[[reduction]] <- AnnotateNMF.DimReduc(object@reductions[[reduction]], object@meta.data)
  }
  object
}

#' Plot metadata enrichment in NMF factors
#' 
#' After running \code{AnnotateNMF}, this function returns a dot plot of the results
#' 
#' @inheritParams AnnotateNMF.DimReduc
#' @param plot.field name of field in \code{meta.data} to plot, if \code{NULL}, plots the first field in the annotation results
#' @return ggplot2 object
#' @aliases AnnotationPlot
#' @rdname AnnotateNMF
#' @export
#' @importFrom stats reshape t.test
#' 
AnnotationPlot.DimReduc <- function(object, plot.field = NULL, ...){
  if(!("annotations" %in% names(object@misc))){
    stop("the ", reduction, " reduction of this object has no 'annotations' slot. Run 'AnnotateNMF' first.")
  }
  if(is.null(field)){
    field <- names(object@misc$annotations)[[1]]
  } else {
    if(!(field %in% names(object@misc$annotations))){
      stop("specified field was not in the annotation object")
    }
    if(length(field) > 1) field <- field[[1]]
  }
  # construct a matrix of pvalues and fc
  ann <- object@misc$annotations[[field]]
  pvals <- reshape(ann, timevar = "group", idvar = "factor", direction = "wide", drop = "fc")
  fc <- reshape(ann, timevar = "group", idvar = "factor", direction = "wide", drop = "p")
  rownames(pvals) <- pvals[,1]
  rownames(fc) <- fc[,1]
  pvals <- pvals[, -1]
  fc <- fc[, -1]
  pvals <- as.matrix(pvals)
  fc <- as.matrix(fc)
  colnames(fc) <- colnames(pvals) <- sapply(colnames(pvals), function(x) substr(x, 3, nchar(x)))
  fc[fc < 0] <- 0
  idx1 <- hclust(dist(t(fc)), method = "ward.D2")$order
  idx2 <- hclust(dist(fc), method = "ward.D2")$order
  fc <- fc[idx2, idx1]
  pvals <- pvals[idx2, idx1]
  fc[fc == 0] <- NA
  pvals <- -log10(pvals)
  pvals[is.infinite(pvals)] <- 100
  pvals[pvals > 100] <- 100
  df <- cbind(reshape2::melt(fc), as.vector(pvals))
  colnames(df) <- c("factor", "field", "fc", "pval")
  df$factor <- factor(df$factor, levels = unique(df$factor))
  suppressWarnings(print(ggplot(df, aes(factor, field, color = pval, size = fc)) + 
                           geom_point() + 
                           scale_color_viridis_c(option = "B", end = 0.9) + 
                           theme_classic() + 
                           labs(y = field, x = "NMF factor", color = "p-value\n(-log10)", size = "fold enrichment") + 
                           theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))))
}

#' @inheritParams AnnotateNMF.Seurat
#' @aliases AnnotationPlot
#' @rdname AnnotateNMF
#' @export
#' 
AnnotationPlot.Seurat <- function(object, plot.field = NULL, reduction = "nmf", ...){
  AnnotationPlot(object@reductions[[reduction]])
}

#' @export
#' @rdname AnnotateNMF
#'
AnnotateNMF <- function(object, ...) {
  UseMethod("AnnotateNMF")
}

#' @export
#' @rdname AnnotateNMF
#'
AnnotationPlot <- function(object, ...) {
  UseMethod("AnnotationPlot")
}

#' @export
#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
.S3method("AnnotationPlot", "Seurat", AnnotationPlot.Seurat)


#' @export
#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
.S3method("AnnotateNMF", "Seurat", AnnotateNMF.Seurat)

#' @export
#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
.S3method("AnnotationPlot", "DimReduc", AnnotationPlot.DimReduc)


#' @export
#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
.S3method("AnnotateNMF", "DimReduc", AnnotateNMF.DimReduc)


#' @export
#' @rdname RankPlot
#' @param reduction the NMF reduction slot name (result of \code{RunNMF} where \code{k} was an array)
#' @param detail.level of detail to plot, \code{1} for test set reconstruction error at convergence of each factorization, \code{2} for test set reconstruction error at each fitting iteration of each factorization
RankPlot.Seurat <- function(object, reduction = "nmf", detail.level = 1, ...) {
  if (detail.level == 2) {
    plot(subset(object@reductions[[reduction]]@misc$cv_data, iter >= 5), detail.level)
  } else {
    plot(object@reductions[[reduction]]@misc$cv_data, detail.level)
  }
}

#' Plot NMF cross-validation results
#'
#' Given a NMF reduction at multiple ranks, plot rank vs. test set reconstruction error to determine the optimal rank.
#'
#' @param object a Seurat object or a \code{data.frame} that is the result of \code{RunNMF}
#' @param reduction name of the NMF reduction in the Seurat object (result of \code{RunNMF}) for which multiple \code{ranks} were computed.
#' @param ... not implemented
#' @export
#' @return A ggplot2 object
#' @aliases RankPlot
#'
RankPlot <- function(object, reduction = "nmf", ...) {
  UseMethod("RankPlot")
}

#' @export
#' @rdname RunNMF
#'
RunNMF <- function(object, ...) {
  UseMethod("RunNMF")
}

#' @export
#' @rdname ProjectData
#'
ProjectData <- function(object, ...) {
  UseMethod("ProjectData")
}

#' @export
#' @rdname RunLNMF
#'
RunLNMF <- function(object, ...) {
  UseMethod("RunLNMF")
}

#' @export
#' @rdname RunLNMF
#'
MetadataPlot <- function(object, ...) {
  UseMethod("MetadataPlot")
}

#' @export
#' @rdname RunNMF
#' @name RunNMF
#'
.S3method("RunNMF", "Seurat", RunNMF.Seurat)

#' Run LNMF on a Seurat object
#'
#' S3 method for Seurat that runs the \code{singlet::RunLNMF} function.
#'
#' @method RunLNMF Seurat
#' @export
#' @rdname RunLNMF
#' @name RunLNMF
#'
.S3method("RunLNMF", "Seurat", RunLNMF.Seurat)

#' Plot NMF cross-validation results given a Seurat object
#'
#' S3 method for Seurat that runs the \code{singlet::RunNMF} function.
#'
#' @method RankPlot Seurat
#' @export
#' @rdname RankPlot
#' @name RankPlot
#'
.S3method("RankPlot", "Seurat", RankPlot.Seurat)

#' @export
#' @rdname RunLNMF
#' @name MetadataPlot
#'
MetadataPlot.Seurat <- function(object, split.by, reduction = "lnmf", ...) {
  if (!(reduction %in% names(object@reductions))) {
    stop("this Seurat object does not contain the requested reductions slot")
  }
  plot(MetadataSummary(t(object@reductions[[reduction]]@cell.embeddings), object@meta.data[[split.by]]))
}

.S3method("MetadataPlot", "Seurat", MetadataPlot.Seurat)

#' @export
#' @rdname RunLNMF
#'
GetSharedFactors <- function(object, split.by, reduction = "lnmf") {
  if (!(reduction %in% names(object@reductions))) {
    stop("this Seurat object does not contain the requested reductions slot")
  }
  # which(rowSums(object@reductions[[reduction]]@misc$link_matrix == 0) == 0)
  which(!(colnames(object@reductions[[reduction]]@cell.embeddings) %in% names(which(apply(MetadataSummary(t(object@reductions[[reduction]]@cell.embeddings), object@meta.data[[split.by]]), 2, function(x) min(x) == 0)))))
}

#' @export
#' @rdname RunLNMF
#'
GetUniqueFactors <- function(object, split.by, reduction = "lnmf") {
  if (!(reduction %in% names(object@reductions))) {
    stop("this Seurat object does not contain the requested reductions slot")
  }
  # which(rowSums(object@reductions[[reduction]]@misc$link_matrix == 0) > 0)
  which((colnames(object@reductions[[reduction]]@cell.embeddings) %in% names(which(apply(MetadataSummary(t(object@reductions[[reduction]]@cell.embeddings), object@meta.data[[split.by]]), 2, function(x) min(x) == 0)))))
}

#' Run Gene Set Enrichment Analysis on a Reduction
#'
#' Run GSEA to identify gene sets that are enriched within NMF factors.
#'
#' @import msigdbr fgsea
#' @param object a Seurat object
#' @param reduction dimensional reduction to use
#' @param species species for which to load gene sets
#' @param category msigdbr gene set category (i.e. "H", "C5", etc.)
#' @param min.size minimum number of terms in a gene set
#' @param max.size maximum number of terms in a gene set
#' @param dims factors in the reduction to use, default \code{NULL} for all factors
#' @param verbose print progress to console
#' @param padj.sig significance cutoff for BH-adjusted p-values (default 0.01)
#' @returns a Seurat object, with GSEA information in the misc slot. BH-adj p-values are on a -log10 scale.
#' @export
#'
RunGSEA <- function(object, reduction = "nmf", species = "Homo sapiens", category = "C5",
                    min.size = 10, max.size = 500, dims = NULL,
                    verbose = TRUE, padj.sig = 0.01) {
  if (verbose) cat("fetching gene sets\n")
  gene_sets <- msigdbr(species = species, category = category)

  if (verbose) cat("filtering pathways\n")
  pathways <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
  pathways <- pathways[lapply(pathways, length) > min.size]

  if (verbose) cat("filtering genes in pathways to those in reduction\n")
  genes_in_pathways <- unique(unlist(pathways))
  w <- object@reductions[[reduction]]@feature.loadings
  if (!is.null(dims)) {
    w <- w[, dims]
  }
  if (verbose) cat("filtering genes in reduction to those in pathways\n")
  w <- w[which(rownames(w) %in% genes_in_pathways), ]
  pathways <- lapply(pathways, function(x) x[x %in% rownames(w)])
  v <- lapply(pathways, length)
  pathways <- pathways[which(v > min.size & v < max.size)]

  cat("running GSEA on", ncol(w), "factors...\n")
  pb <- utils::txtProgressBar(min = 0, max = ncol(w), style = 3)
  results <- list()
  for (i in 1:ncol(w)) {
    ranks <- sort(w[, i])
    results[[i]] <- suppressWarnings(fgseaMultilevel(
      pathways, ranks,
      minSize = min.size, maxSize = max.size, scoreType = "pos"
    ))
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  pval <- do.call(cbind, lapply(results, function(x) x$pval))
  padj <- do.call(cbind, lapply(results, function(x) x$padj))
  es <- do.call(cbind, lapply(results, function(x) x$ES))
  nes <- do.call(cbind, lapply(results, function(x) x$NES))
  rownames(pval) <- rownames(padj) <- rownames(es) <- rownames(nes) <- results[[1]]$pathway

  idx <- which(apply(padj, 1, function(x) min(x) < padj.sig))

  if (!is.null(dims)) {
    dims <- paste0(reduction, dims)
  } else {
    dims <- paste0(reduction, 1:ncol(object@reductions[[reduction]]))
  }
  colnames(pval) <- colnames(padj) <- colnames(es) <- colnames(nes) <- dims

  # reorder with hclust
  padj <- -log10(padj)
  pval <- -log10(pval)
  row_order <- hclust(dist(padj), method = "ward.D2")$order
  col_order <- hclust(dist(t(padj)), method = "ward.D2")$order
  pval <- pval[row_order, col_order]
  padj <- padj[row_order, col_order]
  es <- es[row_order, col_order]
  nes <- nes[row_order, col_order]
  object@reductions[[reduction]]@misc$gsea <- list("pval" = pval, "padj" = padj, "es" = es, "nes" = nes)
  object
}

#' Plot GSEA results on a heatmap
#'
#' Plot top GSEA terms for each NMF factor on a heatmap
#'
#' @param object Seurat object
#' @param reduction a dimensional reduction for which GSEA analysis has been performed
#' @param max.terms.per.factor show this number of top terms for each factor
#' @return ggplot2 object
#' @export
GSEAHeatmap <- function(object, reduction = "nmf", max.terms.per.factor = 3) {
  df <- object@reductions[[reduction]]@misc$gsea$padj

  # find markers for each factor based on the proportion of signal in that factor
  df2 <- as.matrix(Diagonal(x = 1 / rowSums(df)) %*% df)
  terms <- c()
  for (i in 1:ncol(df2)) {
    terms_i <- df[, i]
    idx <- terms_i > -log10(0.05)
    terms_i <- terms_i[idx]
    terms_j <- df2[idx, i]
    v <- sort(terms_j, decreasing = TRUE)
    if (length(v) > max.terms.per.factor) {
      terms <- c(terms, names(v)[1:max.terms.per.factor])
    } else {
      terms <- c(terms, names(v))
    }
  }
  terms <- unique(terms)
  df <- df[terms, ]

  rownames(df) <- sapply(rownames(df), function(x) {
    ifelse(nchar(x) > 48, paste0(substr(x, 1, 45), "..."), x)
  })

  # remove terms that are broadly significant
  v <- which((rowSums(df > -log10(0.05)) > (ncol(df) / 2)))
  if (length(v) > 0) {
    df <- df[-v, ]
  }
  df <- reshape2::melt(df)
  ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "B") +
    theme_classic() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      x = "NMF factor",
      y = "GO Term",
      fill = "BH-adj p-value\n(-log10)"
    ) +
    theme(
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
}

