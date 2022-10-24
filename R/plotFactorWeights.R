#' convenience function to map one or more factors along a genome using igvR
#'
#' @param   object    an nmf object or something with a @w weights matrix
#' @param   gr        a GRanges object with coordinates for the features 
#' @param   factors   which factors to plot weights for (default: 1, 2, 3)
#' @param   plot      use igvR to plot the factors? (TRUE, if igvR detected)
#' 
#' @return            the GRanges gr, but with factor weights added as mcols
#'
#' @details
#'  This function presumes a GRanges object will be supplied, which in turn
#'  presumes that the GenomicRanges package is installed from Bioconductor. 
#'  Further, if plot == TRUE, the igvR package is presumed to be installed. 
#'  If either of these presumptions are false, or if factor weights cannot
#'  be mapped to identifiers in the GRanges, this function will fail. 
#' 
#' @export
#'
plotFactorWeights <- function(object, gr, factors=1:3, plot=FALSE) {

  requireNamespace("GenomicRanges")
  stopifnot(is(gr, "GRanges"))
  stopifnot(all(rownames(object@w) %in% names(gr)))
  gr <- gr[rownames(object@w)]
  
  for (fact in factors) {
    if (is.numeric(fact) | is.integer(fact)) fact <- colnames(object@w)[fact]
    mcols(gr)[, fact] <- object@w[, fact]
  }

  if (plot) {
    requireNamespace("igvR")
    message("igvR support is in process")
  }

  return(gr)

}
