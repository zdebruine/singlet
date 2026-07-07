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
#' @param scoreType scoring type passed to \code{fgsea::fgsea}: \code{"pos"} (default, positive weights only), \code{"std"} or \code{"neg"}
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
RunGSEA <- function(object,...){
  UseMethod("RunGSEA")
}

#' Annotate NMF model with cell or sample metadata
#'
#' @rdname RunGSEA
#' @aliases RunGSEA
#'
#' @export
#'
RunGSEA.Seurat = function(object, ID = "gene_symbol", reduction = "nmf", species = "Homo sapiens", category = "C5",
                          min.size = 10, max.size = 500, dims = NULL,
                          verbose = TRUE, padj.sig = 0.01, add.noise = FALSE, scoreType = "pos", gsea.name = "gsea",remove.zeros.per.factor = FALSE,...){
  
  w <- object@reductions[[reduction]]@feature.loadings
  
  object@reductions[[reduction]]@misc[[gsea.name]] = RunGSEA(object = w,
                                                                    ID = ID,
                                                                    species = species,
                                                                    category = category,
                                                                    min.size = min.size,
                                                                    max.size = max.size,
                                                                    dims = dims,
                                                                    verbose = verbose,
                                                                    padj.sig = padj.sig,
                                                                    add.noise = add.noise,
                                                                    scoreType = scoreType,
                                                                    remove.zeros.per.factor = remove.zeros.per.factor,
                                                                    ...)
  
  return(object)
}

#' @rdname RunGSEA
#' @name RunGSEA
#'
#' @export
#'
.S3method("RunGSEA", "Seurat", RunGSEA.Seurat)


#' Annotate NMF model with cell or sample metadata
#'
#' @rdname RunGSEA
#' @aliases RunGSEA
#'
#' @export
#'
RunGSEA.SingleCellExperiment = function(object, ID = "gene_symbol", reduction = "nmf", species = "Homo sapiens", category = "C5",
                                        min.size = 10, max.size = 500, dims = NULL,
                                        verbose = TRUE, padj.sig = 0.01, add.noise = FALSE, scoreType = "pos", gsea.name = "gsea",remove.zeros.per.factor = FALSE,...){
  object@metadata[[reduction]] = RunGSEA(object = object@metadata[[reduction]],
                                         ID = ID,
                                         species = species,
                                         category = category,
                                         min.size = min.size,
                                         max.size = max.size,
                                         dims = dims,
                                         verbose = verbose,
                                         padj.sig = padj.sig,
                                         add.noise = add.noise,
                                         scoreType = scoreType,
                                         gsea.name = gsea.name,
                                         remove.zeros.per.factor = remove.zeros.per.factor,
                                         ...)
  
  return(object)
}

#' @rdname RunGSEA
#' @name RunGSEA
#'
#' @export
#'
.S3method("RunGSEA", "SingleCellExperiment", RunGSEA.SingleCellExperiment)

#' Annotate NMF model with cell or sample metadata
#'
#' @rdname RunGSEA
#' @aliases RunGSEA
#'
#' @export
#'
RunGSEA.nmfList = function(object, ID = "gene_symbol", species = "Homo sapiens", category = "C5",
                        min.size = 10, max.size = 500, dims = NULL,
                        verbose = TRUE, padj.sig = 0.01, add.noise = FALSE, scoreType = "pos", gsea.name = "gsea",remove.zeros.per.factor = FALSE,...){
  w = object$w
  
  object[[gsea.name]] = RunGSEA(object = w,
                                       ID = ID,
                                       species = species,
                                       category = category,
                                       min.size = min.size,
                                       max.size = max.size,
                                       dims = dims,
                                       verbose = verbose,
                                       padj.sig = padj.sig,
                                       add.noise = add.noise,
                                       scoreType = scoreType,
                                       remove.zeros.per.factor = remove.zeros.per.factor,
                                       ...)
  return(object)
}

#' @rdname RunGSEA
#' @name RunGSEA
#'
#' @export
#'
.S3method("RunGSEA", "nmfList", RunGSEA.nmfList)

#' Annotate NMF model with cell or sample metadata
#'
#' @rdname RunGSEA
#' @aliases RunGSEA
#'
#' @export
#'
RunGSEA.nmf = function(object, ID = "gene_symbol", species = "Homo sapiens", category = "C5",
                       min.size = 10, max.size = 500, dims = NULL,
                       verbose = TRUE, padj.sig = 0.01, add.noise = FALSE, scoreType = "pos", gsea.name = "gsea",remove.zeros.per.factor = FALSE,...){
  w <- object@w
  
  object@misc[[gsea.name]] = RunGSEA(object = w,
                                            ID = ID,
                                            species = species,
                                            category = category,
                                            min.size = min.size,
                                            max.size = max.size,
                                            dims = dims,
                                            verbose = verbose,
                                            padj.sig = padj.sig,
                                            add.noise = add.noise,
                                            scoreType = scoreType,
                                            remove.zeros.per.factor = remove.zeros.per.factor,
                                            ...)
  return(object)
}

#' @rdname RunGSEA
#' @name RunGSEA
#'
#' @export
#'
.S3method("RunGSEA", "nmf", RunGSEA.nmf)



#' Annotate NMF model with cell or sample metadata
#'
#' @rdname RunGSEA
#' @aliases RunGSEA
#'
#' @export
#'
RunGSEA.matrix <- function(object, ID = "gene_symbol", species = "Homo sapiens", category = "C5",
                           min.size = 10, max.size = 500, dims = NULL,
                           verbose = TRUE, padj.sig = 0.01, add.noise = FALSE, scoreType = "pos",remove.zeros.per.factor = FALSE,...) {
  
  if (verbose) cat("fetching gene sets\n")
  gene_sets <- msigdbr(species = species, category = category, ...)
  
  if (verbose) cat("filtering pathways\n")
  pathways <- split(x = gene_sets[[ID]], f = gene_sets$gs_name)
  pathways <- pathways[lapply(pathways, length) > min.size]
  
  if (verbose) cat("filtering genes in pathways to those in reduction\n")
  genes_in_pathways <- unique(unlist(pathways))
  
  if (!is.null(dims)) object <- object[, dims,drop=FALSE]
  
  
  #this might be redundant fgsea looks to already be removing pathways in preparePathways() based on genes
  #they are also removing pathways that do not meet the max and min size
  #https://github.com/alserglab/fgsea/blob/d6f1670c60bbc48ad1a95f98dc52921ec98d2542/R/util.R#L1
  if (verbose) cat("filtering genes in reduction to those in pathways\n")
  object <- object[which(rownames(object) %in% genes_in_pathways), ,drop=FALSE]
  pathways <- lapply(pathways, function(x) x[x %in% rownames(object)])
  v <- lapply(pathways, length)
  pathways <- pathways[which(v > min.size & v < max.size)]
  
  #remove genes that are 0 across all factors
  #abs make it work for more than just nmf factors
  row_sum = rowSums(abs(object))
  if (any(row_sum == 0)){
    cat("Rows with all 0 detected. Removing...\n")
    object = object[row_sum !=0,,drop=FALSE]
  }
  
  
  cat("running GSEA on", ncol(object), "factors...\n")
  pb <- utils::txtProgressBar(min = 0, max = ncol(object), style = 3)
  results <- list()
  ranksList = list()
  for (i in 1:ncol(object)) {
    #this is also redundant as fgsea also sorts automatically in preparePathwaysAndStats()
    #https://github.com/alserglab/fgsea/blob/d6f1670c60bbc48ad1a95f98dc52921ec98d2542/R/fgsea.R#L41
    ranks <- sort(object[, i])
    
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
    ranksList[[i]] = ranks
    
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  factor_names = colnames(object)
  if (is.null(factor_names)) {
    factor_names = paste0("factor", seq_len(ncol(object)))
  }
  names(results) = factor_names
  names(ranksList) = factor_names
  
  
  return(list(gsea = results,ranks = ranksList,pathways = pathways))
}
#' @rdname RunGSEA
#' @name RunGSEA
#'
#' @export
#'
.S3method("RunGSEA", "matrix", RunGSEA.matrix)



