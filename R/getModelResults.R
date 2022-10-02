#' extract data.frame of lods and pvalues for differential factor representation
#' 
#' log-odds of non-null differences for a response by a factor are in fit$lods
#' (which will usually be a matrix), and one-sided p-values for the moderated t
#' test are computed from fit$t and fit$df.total using pt(t, df, lower=FALSE),
#' then adjusted using the step-up procedure of Benjamini & Hochberg.
#' 
#' @param fit   an lmFit result from limma, shrunken with eBayes()
#'
#' @return      a data.frame with columns 'factor', 'group', 'fc', and 'p' 
#' 
#' @importFrom  reshape2 melt
#' @import      limma
#'
#' @export
getModelResults <- function(fit) { 

  # fits are centered, so use signed lods for evidence 
  fc <- reshape2::melt(fit$lods * sign(fit$coefficients))
  names(fc) <- c("factor", "group", "fc")

  # computed BH-adjusted moderated one-way p-values 
  p <- pt(fit$t, fit$df.total, lower=FALSE)
  padj <- matrix(p.adjust(p, method="fdr"), dimnames=dimnames(p), ncol=ncol(p))
  fdr <- reshape2::melt(padj)
  names(fdr) <- c("factor", "group", "p")

  # now bolt them back together
  merge(fc, fdr)

}
