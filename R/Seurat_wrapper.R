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
#' @param test.set.density approximate density of the test set (default 0.05)
#' @param tol tolerance of the fit (correlation distance of the model across consecutive iterations)
#' @param maxit maximum number of fitting iterations
#' @param L1 L1/LASSO penalty to increase sparsity of the model
#' @param threads number of threads to use (0 = let OpenMP use all available threads)
#' @param split.by Attribute in \code{@metadata} for splitting, if applicable. Data will be weighted such that each group contributes equally to the NMF model.
#' @param precision \code{"double"} or \code{"float"} for numerical precision. \code{"float"} may be faster, but numerical instability may result in more iterations to achieve desired tolerances.
#' @details Use \code{set.seed()} to guarantee reproducibility!
#' @export
#' @aliases RunNMF
#' @seealso \code{\link{RunLNMF}}, \code{\link{RankPlot}}, \code{\link{MetadataSummary}}
#' @return Returns a Seurat object with the NMF model stored in the reductions slot
#' 
RunNMF.Seurat <- function(
  object, 
  assay = NULL, 
  reduction.name = "nmf", 
  reduction.key = "NMF_", 
  verbose = 2, 
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
  v <- object@assays[[assay]]@data@x
  if(sum(as.integer(v)) == sum(v)){
    object <- Seurat::NormalizeData(object, assay = assay)
  }
  A <- object@assays[[assay]]@data
  rnames <- rownames(A)
  cnames <- colnames(A)
  
  if(!is.null(split.by)){
    split.by <- as.integer(as.numeric(as.factor(object@meta.data[[split.by]]))) - 1
    if(any(is.na(split.by)))
      stop("'split.by' cannot contain NA values")
    A <- weight_by_split(A, split.by, length(unique(split.by)))
  }
  At <- Matrix::t(A)
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
#' @inheritParams RunNMF.Seurat
#' @param split.by attribute in \code{@metadata} for splitting. Data will be weighted such that each group contributes equally to the LNMF model.
#' @param link.cutoff if the relative contribution of samples in any given group to a factor falls below \code{link.cutoff}, unlink it from the factor. \code{link.cutoff = 1} means a factor must contribute exactly equally before being unlinked.
#' @param link.balance after initial linking step, weight all shared factors such that each dataset is represented as equally as possible in each factor.
#' @param link.balance.tol relative change in factor representation within sample groups between balancing iterations at which to call convergence.
#' @param link.balance.rate proportion of difference between current factor weight and equal representation of factors in a sample group to target in a single iteration (default 0.1).  
#' @param balance.maxit maximum number of iterations for the balancing step
#' @param reduction.use NMF reduction to use for initializing the linked factorization.
#' @param reduction.name name to store resulting DimReduc object as
#' @param reduction.key key for resulting DimReduc
#' @param verbose print fitting progress to console
#' @param tol tolerance of the fit (correlation distance of the model across consecutive iterations)
#' @param maxit maximum number of fitting iterations
#' @param L1 L1/LASSO penalty to increase sparsity of the model
#' @param threads number of threads to use (0 = let OpenMP use all available threads)
#' @param reduction reduction to use for metadata analysis
#' @seealso \code{\link{RunNMF}}, \code{\link{RankPlot}}, \code{\link{MetadataSummary}}
#' @export
#' @details Use \code{set.seed()} to guarantee reproducibility!
#' @aliases RunLNMF
#' @rdname RunLNMF
#' @return Returns a Seurat object with the NMF model stored in the reductions slot
#' 
RunLNMF.Seurat <- function(
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
      A <- weight_by_split(t(A), split.by, length(unique(split.by)))
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
    m <- MetadataSummary(lnmf_model$h, split.by, FALSE)
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
      m <- MetadataSummary(lnmf_model$h, split.by, FALSE)
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

#' @export
#' @rdname RankPlot
#' @param reduction the NMF reduction slot name (result of \code{RunNMF} where \code{k} was an array)
RankPlot.Seurat <- function(object, reduction = "nmf"){
  plot(object@reductions[[reduction]]@misc$cv_data)
}

#' Plot NMF cross-validation results
#' 
#' Given a NMF reduction at multiple ranks, plot rank vs. test set reconstruction error to determine the optimal rank.
#'
#' @param object a Seurat object or a \code{data.frame} that is the result of \code{RunNMF}
#' @param reduction name of the NMF reduction in the Seurat object (result of \code{RunNMF}) for which multiple \code{ranks} were computed.
#'
#' @export
#' @return A ggplot2 object
#' @aliases RankPlot
#' 
RankPlot <- function(object, reduction = "nmf") {
  UseMethod("RankPlot")
}

#' @export 
#' @rdname RunNMF
#'
RunNMF <- function(object, ...){
  UseMethod("RunNMF")
}

#' @export
#' @rdname RunLNMF
#'
RunLNMF <- function(object, ...){
  UseMethod("RunLNMF")
}

#' @export
#' @rdname RunLNMF
#'
MetadataPlot <- function(object, ...){
  UseMethod("MetadataPlot")
}

#' @export
#' @rdname RunNMF
#'
.S3method("RunNMF", "Seurat", RunNMF.Seurat)

#' @export
#' @rdname RunLNMF
#'
.S3method("RunLNMF", "Seurat", RunLNMF.Seurat)

#' @export
#' @rdname RankPlot
#'
.S3method("RankPlot", "Seurat", RankPlot.Seurat)

#' @export
#' @rdname RunLNMF
#'
MetadataPlot.Seurat <- function(object, split.by, reduction = "lnmf"){
  if(!(reduction %in% names(object@reductions)))
    stop("this Seurat object does not contain the requested reductions slot")
  plot(MetadataSummary(t(object@reductions[[reduction]]@cell.embeddings), object@meta.data[[split.by]]))
}

.S3method("MetadataPlot", "Seurat", MetadataPlot.Seurat)

#' @export
#' @rdname RunLNMF
#'
GetSharedFactors <- function(object, split.by, reduction = "lnmf"){
  if(!(reduction %in% names(object@reductions)))
    stop("this Seurat object does not contain the requested reductions slot")
  which(!(colnames(object@reductions[[reduction]]@cell.embeddings) %in% names(which(apply(MetadataSummary(t(object@reductions[[reduction]]@cell.embeddings), object@meta.data[[split.by]]), 1, function(x) min(x) == 0)))))
}

#' @export
#' @rdname RunLNMF
#'
GetUniqueFactors <- function(object, split.by, reduction = "lnmf"){
  if(!(reduction %in% names(object@reductions)))
    stop("this Seurat object does not contain the requested reductions slot")
  which((colnames(object@reductions[[reduction]]@cell.embeddings) %in% names(which(apply(MetadataSummary(t(object@reductions[[reduction]]@cell.embeddings), object@meta.data[[split.by]]), 1, function(x) min(x) == 0)))))
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
#' @param collapse filter pathways to remove highly redundant terms
#' @param dims factors in the reduction to use, default \code{NULL} for all factors
#' @param verbose print progress to console
#' @param padj.size significance cutoff for BH-adjusted p-values (default 0.01)
#' @param ... additional arguments to \code{fgseaMultilevel}
#' @returns a Seurat object, with GSEA information in the misc slot. BH-adj p-values are on a -log10 scale.
#' @export
#' 
RunGSEA <- function(object, reduction = "nmf", species = "Homo sapiens", category = "C5", 
                    min.size = 10, max.size = 500, collapse = TRUE, dims = NULL, 
                    verbose = TRUE, padj.sig = 0.01, ...){
  
  if(verbose) cat("fetching gene sets\n")
  gene_sets = msigdbr(species = species, category = category)
  
  if(verbose) cat("filtering pathways\n")
  pathways = split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
  pathways <- pathways[lapply(pathways, length) > min.size]
  
  if(verbose) cat("filtering genes in pathways to those in reduction\n")
  genes_in_pathways <- unique(unlist(pathways))
  w <- object@reductions[[reduction]]@feature.loadings
  if(!is.null(dims))
    w <- w[,dims]
  if(verbose) cat("filtering genes in reduction to those in pathways\n")
  w <- w[which(rownames(w) %in% genes_in_pathways), ]
  pathways <- lapply(pathways, function(x) x[x %in% rownames(w)])
  v <- lapply(pathways, length)
  pathways <- pathways[which(v > min.size & v < max.size)]
  
  cat("running GSEA on", ncol(w), "factors...\n")
  pb <- utils::txtProgressBar(min = 0, max = ncol(w), style = 3)
  results <- collapsed <- list()
  for(i in 1:ncol(w)){
    ranks <- sort(w[, i])
    results[[i]] <- suppressWarnings(fgseaMultilevel(
      pathways, ranks, minSize = min.size, maxSize = max.size, scoreType = "pos", ...))
    
    if(collapse){
      collapsedPathways <- collapsePathways(
        results[[i]][order(pval)][padj < padj.sig], pathways, ranks)
      collapsed[[i]] <- results[[i]][pathway %in% collapsedPathways$mainPathways]
    }
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  # filter results to only collapsed pathways
  if(collapse){
    collapsed <- unique(unlist(collapsed))
    results <- lapply(results, function(x) {
      x[pathway %in% collapsed, ]
    })
  }
  
  pval <- do.call(cbind, lapply(results, function(x) x$pval))
  padj <- do.call(cbind, lapply(results, function(x) x$padj))
  es <- do.call(cbind, lapply(results, function(x) x$ES))
  nes <- do.call(cbind, lapply(results, function(x) x$NES))
  rownames(pval) <- rownames(padj) <- rownames(es) <- rownames(nes) <- results[[1]]$pathway

  idx <- which(apply(padj, 1, function(x) min(x) < padj.sig))
  
  if(!is.null(dims)) {
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
GSEAHeatmap <- function(object, reduction = "nmf", max.terms.per.factor = 3){
  df <- object@reductions[[reduction]]@misc$gsea$padj
  
  # find markers for each factor based on the proportion of signal in that factor
  df2 <- as.matrix(Diagonal(x = 1 / rowSums(df)) %*% df)
  terms <- c()
  for(i in 1:ncol(df2)){
    terms_i <- df[,i]
    idx <- terms_i > -log10(0.05)
    terms_i <- terms_i[idx]
    terms_j <- df2[idx, i]
    v <- sort(terms_j, decreasing = TRUE)
    if(length(v) > max.terms.per.factor){
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
  if(length(v) > 0){
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
      fill = "BH-adj p-value\n(-log10)") + 
    theme(
      axis.text.y = element_text(size = 6), 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

ClusteringToDimReduc <- function(object, group.by = "ident"){

  ifelse(ident == "active.ident", idents <- as.vector(Idents(object)), idents <- as.vector(object@meta.data[[ident]]))
  ident.names <- unique(idents)
  if (verbose > 0) pb <- txtProgressBar(char = "=", style = 3, max = length(ident.names), width = 50)
  m <- list()
  for (i in 1:length(ident.names)) {
    m[[i]] <- Matrix::rowMeans(data[, which(idents == ident.names[i])])
    if (verbose > 0) setTxtProgressBar(pb = pb, value = i)
  }
  result <- do.call(cbind, m)
  colnames(result) <- ident.names
}