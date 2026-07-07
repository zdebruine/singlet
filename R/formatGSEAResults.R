#' Format GSEA Results into Matrices
#'
#' Combine GSEA results across factors into aligned matrices of enrichment
#' statistics. Optionally filter to a specified set of pathways.
#'
#' @param gsea A result object returned by runGSEA().
#' @param pathways Optional character vector of pathway names to retain.
#'   Defaults to NULL, retaining all pathways.
#' @param log10pval Logical indicating whether to apply negative log10
#'   transformation to p-values. Default TRUE.
#'
#' @return A list containing four matrices: \code{pval}, \code{padj},
#'   \code{es}, and \code{nes}. Rows are pathways and columns are factors.
#'   Missing pathway results are represented by NA.
#'
#' @keywords internal
formatGSEA = function(gsea, pathways = NULL, log10pval = TRUE) {
  
  results = gsea$gsea
  factor_names = names(results)
  
  if (is.null(results) || length(results) == 0) {
    stop("gsea must contain a non-empty list of GSEA results")
  }
  
  if (is.null(factor_names)) {
    stop("GSEA results must have factor names")
  }
  
  if (!is.null(pathways)) {
    results = lapply(results, function(x) {
      x[x$pathway %in% pathways, , drop = FALSE]
    })
  }
  
  if (all(vapply(results, nrow, integer(1)) == 0)) {
    
    empty = matrix(
      numeric(0),
      nrow = 0,
      ncol = length(results),
      dimnames = list(character(0), factor_names)
    )
    
    return(list(
      pval = empty,
      padj = empty,
      es = empty,
      nes = empty
    ))
  }
  
  results = padData(results)
  
  pval = do.call(cbind, lapply(results, `[[`, "pval"))
  padj = do.call(cbind, lapply(results, `[[`, "padj"))
  es = do.call(cbind, lapply(results, `[[`, "ES"))
  nes = do.call(cbind, lapply(results, `[[`, "NES"))
  
  pathway_names = results[[1]]$pathway
  
  rownames(pval) = rownames(padj) =
    rownames(es) = rownames(nes) = pathway_names
  
  colnames(pval) = colnames(padj) =
    colnames(es) = colnames(nes) = factor_names
  
  if (log10pval) {
    pval = -log10(pval)
    padj = -log10(padj)
  }
  
  return(list(
    pval = pval,
    padj = padj,
    es = es,
    nes = nes
  ))
}

#' Collapse Redundant GSEA Pathways
#'
#' Identify representative pathways for each factor using
#' fgsea::collapsePathways().
#'
#' @param gsea A result object returned by runGSEA(), containing
#'   \code{gsea}, \code{ranks}, and \code{pathways}.
#' @param padj.sig Adjusted p-value threshold for selecting pathways to collapse.
#' @param pval.threshold p-value threshold for collapsing pathways
#' @param ... Additional arguments passed to fgsea::collapsePathways().
#'
#' @return A named list of collapsePathways() results, one per factor.
#'   Factors without significant pathways contain NULL.
#' @importFrom fgsea collapsePathways
#' @export
collapseGSEA = function(gsea, padj.sig = 0.05, pval.threshold = 0.05, ...) {
  
  start_time = Sys.time()
  message("Collapsing GSEA terms...")
  
  results = gsea$gsea
  ranks = gsea$ranks
  pathways = gsea$pathways
  
  if (is.null(results) || is.null(ranks) || is.null(pathways)) {
    stop("gsea must contain 'gsea', 'ranks', and 'pathways'")
  }
  
  if (length(results) == 0 || length(results) != length(ranks)) {
    stop("GSEA results and ranks must have the same nonzero length")
  }
  
  if (!is.numeric(padj.sig) || length(padj.sig) != 1 ||
      is.na(padj.sig) || padj.sig < 0 || padj.sig > 1) {
    stop("padj.sig must be between 0 and 1")
  }
  
  collapsed = lapply(seq_along(results), function(i) {
    
    res = results[[i]]
    
    sig = res[
      !is.na(res$padj) & res$padj < padj.sig,
    ]
    
    if (nrow(sig) == 0) {
      return(NULL)
    }
    
    cc = fgsea::collapsePathways(
      fgseaRes = sig,
      pathways = pathways,
      stats = ranks[[i]],
      pval.threshold = pval.threshold,
      ...
    )
    cc[["padj.sig"]] = padj.sig
    cc[["pval.threshold"]] = pval.threshold
    
    return(cc)
  })
  
  names(collapsed) = names(results)
  
  elapsed = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  message(sprintf("Done (%.1f seconds)", elapsed))
  return(collapsed)
}

#' Extract Representative Pathways
#'
#' Extract the union of representative pathways across factors from
#' collapseGSEA() results.
#'
#' @param collapsed A named list returned by collapseGSEA().
#'
#' @return A character vector of unique representative pathway names.
#'
#' @keywords internal
getMainPathways = function(collapsed) {
  
  main_pathways = lapply(collapsed, function(x) {
    if (is.null(x)) return(character(0))
    x$mainPathways
  })
  
  return(unique(unlist(main_pathways, use.names = FALSE)))
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