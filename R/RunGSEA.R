#' Run Gene Set Enrichment Analysis on a Reduction
#'
#' Run GSEA to identify gene sets that are enriched within NMF factors.
#'
#' @param object a Seurat object
#' @param reduction dimensional reduction to use
#' @param species species for which to load gene sets
#' @param category msigdbr gene set category (i.e. "H", "C5", etc.)
#' @param min.size minimum number of terms in a gene set
#' @param max.size maximum number of terms in a gene set
#' @param dims factors in the reduction to use, default \code{NULL} for all factors
#' @param verbose print progress to console
#' @param padj.sig significance cutoff for BH-adjusted p-values (default 0.01)
#' @return a Seurat object, with GSEA information in the misc slot. BH-adj p-values are on a -log10 scale.
#'
#' @importFrom fgsea msigdbr
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
