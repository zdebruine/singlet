#' convenience function to plot one or more factors along a genome using igvR
#'
#' @param   object    an nmf object or something with a @w weights matrix
#' @param   gr        a GRanges object with coordinates for the features 
#' @param   factors   which factors to plot weights for (default: 1, 2, 3)
#' @param   plot      use igvR to plot the factors? (TRUE, if igvR detected)
#' 
#' @return            the GRanges gr, but with factor weights added as mcols
#'
#' @export
#'
plotFactorWeights <- function(object, gr, factors=1:3, plot=NULL) {
  requireNamespace("GenomicRanges")
  stopifnot(is(gr, "GRanges"))
  stopifnot(all(rownames(object@w) %in% names(gr)))
  gr <- gr[rownames(object@w)]
  
  for (fact in factors) {
    if (is.numeric(fact) | is.integer(fact)) fact <- colnames(object@w)[fact]
    mcols(gr)[, fact] <- object@w[, fact]
  }

  # use :: not require()
  #  if (is.null(plot)) plot <- require(igvR)
  #  if (plot) message("igvR support is in process")

  return(gr)

}
