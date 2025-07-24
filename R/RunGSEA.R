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
#' @param gsea.name gsea name, gsea by default
#' @param remove.zeros.per.factor option to remove genes that have weight 0 within factor (not recommended unless you can prove this is unbiased filtering)
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
                    verbose = TRUE, padj.sig = 0.01, add.noise = FALSE, scoreType = "pos", gsea.name = "gsea",remove.zeros.per.factor = FALSE,...) {

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
  
  
  #this might be redundant fgsea looks to already be removing pathways in preparePathways() based on genes
  #they are also removing pathways that do not meet the max and min size
  #https://github.com/alserglab/fgsea/blob/d6f1670c60bbc48ad1a95f98dc52921ec98d2542/R/util.R#L1
  if (verbose) cat("filtering genes in reduction to those in pathways\n")
  w <- w[which(rownames(w) %in% genes_in_pathways), ]
  pathways <- lapply(pathways, function(x) x[x %in% rownames(w)])
  v <- lapply(pathways, length)
  pathways <- pathways[which(v > min.size & v < max.size)]

  #remove genes that are 0 across all factors
  row_sum = rowSums(w)
  if (any(row_sum == 0)){
    cat("Rows with all 0 detected. Removing...\n")
    w = w[row_sum !=0,]
  }
  
  
  cat("running GSEA on", ncol(w), "factors...\n")
  pb <- utils::txtProgressBar(min = 0, max = ncol(w), style = 3)
  results <- list()
  for (i in 1:ncol(w)) {
    #this is also redundant as fgsea also sorts automatically in preparePathwaysAndStats()
    #https://github.com/alserglab/fgsea/blob/d6f1670c60bbc48ad1a95f98dc52921ec98d2542/R/fgsea.R#L41
    ranks <- sort(w[, i])
    
    #I need to check if this is unbiased filtering.
    if(remove.zeros.per.factor){
      ranks = ranks[ranks != 0]

    }
    
    
    # because everything gets reranked, I need to look at whether adding random noise is going to change the ranks.
    # ranks may change if there are duplicate values, this is likely to occur if you use nmf factors since you get lots of 0
    # removing 0 because I can't confirm that fgsea is already doing this. 
    if (add.noise){
      ranks = ranks+rnorm(length(ranks), sd=0.001)
    }
    
    results[[i]] <- suppressWarnings(fgseaMultilevel(
      pathways, ranks,
      minSize = min.size, maxSize = max.size, scoreType = scoreType
    ))
    
    
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  results = padData(results)
  
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
  
  row_order <- hclust(dist(na.omit(padj)), method = "ward.D2")$order
  col_order <- hclust(dist(t(na.omit(padj))), method = "ward.D2")$order
  pval <- pval[row_order, col_order]
  padj <- padj[row_order, col_order]
  es <- es[row_order, col_order]
  nes <- nes[row_order, col_order]

  if (is(object, "Seurat")) {
    object@reductions[[reduction]]@misc[[gsea.name]] <- 
      list("pval" = pval, "padj" = padj, "es" = es, "nes" = nes)
  } else if (is(object, "nmf")) { 
    object@misc[[gsea.name]] <- 
      list("pval" = pval, "padj" = padj, "es" = es, "nes" = nes)
  }

  object
}

uniquePathways = function(resultList){
  vec = character()
  for (i in 1:length(resultList) ){
    vec = c(vec,resultList[[i]]$pathway)
  }
  
  return(unique(vec))
}

fill_missing = function(df,pathways){
  temp = data.frame(pathway = pathways)
  joined = merge(temp, df, by = "pathway", all.x = TRUE)
  
  reordered_ind = match(pathways,joined$pathway)
  joined = joined[reordered_ind,]
  rownames(joined) = 1:nrow(joined)
  
  return(joined)
}

padData = function(resultList){
  pathways = uniquePathways(resultList)
  
  padded_results = lapply(resultList,fill_missing,pathways = pathways)
  
  return(padded_results)
  
}