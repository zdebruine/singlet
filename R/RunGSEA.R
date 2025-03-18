#' Run Gene Set Enrichment Analysis on a Reduction
#'
#' Run GSEA to identify gene sets that are enriched within NMF factors.
#'
#' @param object a Seurat or RcppML::nmf object
#' @param ID gene ID format (i.e. "gene_symbol", "ensembl_gene")
#' @param reduction dimensional reduction to use (if Seurat)
#' @param species species for which to load gene sets
#' @param category msigdbr gene set category (i.e. "H", "C5", etc.)
#' @param min.size minimum number of terms in a gene set
#' @param max.size maximum number of terms in a gene set
#' @param dims factors in the reduction to use, default \code{NULL} for all factors
#' @param verbose print progress to console
#' @param padj.sig significance cutoff for BH-adjusted p-values (default 0.01)
#' @param add.noise  add random noise to ranks to prevent hanging, known bug in fgsea(https://github.com/alserglab/fgsea/issues/151#issuecomment-2088857387)
#' @param ...   additional params to pass to msigdbr
#'
#' @return a Seurat or nmf object, with GSEA information in the misc slot. BH-adj p-values are on a -log10 scale.
#'
#' @import fgsea
#' @import msigdbr
#'
#' @export
#'
RunGSEA <- function(object, ID = "gene_symbol", reduction = "nmf", species = "Homo sapiens", category = "C5",
                    min.size = 10, max.size = 500, dims = NULL,
                    verbose = TRUE, padj.sig = 0.01, add.noise = FALSE,...) {

  if (verbose) cat("fetching gene sets\n")
  gene_sets <- msigdbr(species = species, category = category, ...)

  if (verbose) cat("filtering pathways\n")
  pathways <- split(x = gene_sets[[ID]], f = gene_sets$gs_name)
  pathways <- pathways[lapply(pathways, length) > min.size]

  if (verbose) cat("filtering genes in pathways to those in reduction\n")
  genes_in_pathways <- unique(unlist(pathways))

  # work on RcppML nmf objects too: 
  if (is(object, "Seurat")) {
    w <- object@reductions[[reduction]]@feature.loadings
  } else if (is(object, "nmf")) { 
    w <- object@w
  }
  if (!is.null(dims)) w <- w[, dims]
  
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

    if (add.noise){
      results[[i]] <- suppressWarnings(fgseaMultilevel(
        pathways, ranks+rnorm(length(ranks), sd=0.001),
        minSize = min.size, maxSize = max.size, scoreType = "pos"
      ))
    }else{
      results[[i]] <- suppressWarnings(fgseaMultilevel(
        pathways, ranks,
        minSize = min.size, maxSize = max.size, scoreType = "pos"
      ))
    }
    
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
  } else if (is(object, "Seurat")) {
    dims <- paste0(reduction, 1:ncol(object@reductions[[reduction]]))
  } else if (is(object, "nmf")) {
    dims <- paste0("nmf", 1:ncol(w))
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

  if (is(object, "Seurat")) {
    object@reductions[[reduction]]@misc$gsea <- 
      list("pval" = pval, "padj" = padj, "es" = es, "nes" = nes)
  } else if (is(object, "nmf")) { 
    object@misc$gsea <- 
      list("pval" = pval, "padj" = padj, "es" = es, "nes" = nes)
  }

  object
}
