#' given data in rankable order, squash it into 0-1 and project onto a Gaussian.
#'
#' @param x     data that can be ranked somehow (usually numeric)
#'
#' @return      the same data projected onto (-Inf, +Inf)
#' 
#' @seealso     normitAll for matrices or other rectangular shaped things
#'
#' @example
#'
#' par(mfrow=c(1,2))
#'
#' x <- rgamma(n=10000, shape=3, scale=5)
#' plot(density(x), main="Gamma-distributed")
#'
#' y <- normit(x)
#' plot(density(y), main="normit-transformed")
#' 
#' par(mfrow=c(1,1))
#'
#' @export
normit <- function(x) { 

  notNA <- which(!is.na(x))
  n <- length(unique(x[notNA]))
  x[notNA] <- qnorm((rank(x[notNA])-0.5)/n)
  return(x)

}
