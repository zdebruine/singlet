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

  # log-odds of differential representation
  stopifnot("lods" %in% names(fit))
  lods <- reshape2::melt(fit$lods)
  names(lods) <- c("factor", "group", "fc")

  # BH-adjusted p-values (treated as one huge vector)
  stopifnot("p.value" %in% names(fit))
  pmat <- fit$p.value
  padj <- matrix(p.adjust(pmat, method="fdr"), 
                 dimnames=dimnames(pmat),
                 ncol=ncol(pmat))
  fdr <- reshape2::melt(padj)
  names(fdr) <- c("factor", "group", "p")

  # now bolt them back together
  merge(lods, fdr)

} 
