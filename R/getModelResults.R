#' extract data.frame of lods and pvalues for differential factor representation
#' 
#' log-odds of non-null differences for a response by a factor are in fit$lods
#' (which will usually be a matrix) and p-values will be in fit$p.value (these
#' are raw, so we adjust them to yield a B-H fdr). The results are then merged.
#' 
#' @param fit   an lmFit result from limma, perhaps shrunken with eBayes()
#'
#' @return      a data.frame with columns 'factor', 'group', 'fc', and 'p' 
#' 
#' @importFrom  reshape2 melt
#' @import      limma
#'
#' @export
getModelResults <- function(fit) { 

  stopifnot("lods" %in% names(fit))
  lods <- reshape2::melt(fit$lods)
  names(lods) <- c("factor", "group", "fc")

  stopifnot("p.value" %in% names(fit))
  pBH <- reshape2::melt(apply(fit$p.value, 2, p.adjust, method="BH"))
  names(pBH) <- c("factor", "group", "p")

  merge(lods, pBH)

} 
