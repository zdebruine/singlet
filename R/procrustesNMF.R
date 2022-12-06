#' @title Rescaling of cell and gene features
#'
#' @param Object_X is the reference Object 
#' @param Object_Y is the query object, to be modified to have similiar coordinates.
#' @param scale is used to scale the axes of UMAP coordinates of the Object_Y.
#' @param symmetric to confirm the symmetric of the rotation.
#' @export


setGeneric("procrustesNMF", function(Object_X,Object_Y,...) {
  standardGeneric("procrustesNMF")
})

setMethod("procrustesNMF", signature("Seurat"), function(Object_X, Object_Y, scale = TRUE, symmetric = FALSE,...) { 
  X <- Object_X
  Y <- Object_Y
  X <- scores(X, display = scores, ...)
  Y <- scores(Y, display = scores, ...)
  if (nrow(X) != nrow(Y))
    stop(gettextf("matrices have different number of rows: %d and %d",
                  nrow(X), nrow(Y)))
  if (ncol(X) < ncol(Y)) {
    warning("X has fewer axes than Y: X adjusted to comform Y\n")
    addcols <- ncol(Y) - ncol(X)
    for (i in 1:addcols) X <- cbind(X, 0)
  }
  ctrace <- function(MAT) sum(MAT^2)
  c <- 1
  if (symmetric) {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
    X <- X/sqrt(ctrace(X))
    Y <- Y/sqrt(ctrace(Y))
  }
  xmean <- apply(X, 2, mean)
  ymean <- apply(Y, 2, mean)
  if (!symmetric) {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
  }
  XY <- crossprod(X, Y)
  sol <- svd(XY)
  A <- sol$v %*% t(sol$u)
  if (scale) {
    c <- sum(sol$d)/ctrace(Y)
  }
  Yrot <- c * Y %*% A
  ## Translation (b) needs scale (c) although Mardia et al. do not
  ## have this. Reported by Christian Dudel.
  b <- xmean - c * ymean %*% A
  R2 <- ctrace(X) + c * c * ctrace(Y) - 2 * c * sum(sol$d)
  reslt <- list(Yrot = Yrot, X = X, ss = R2, rotation = A,
                translation = b, scale = c, xmean = xmean,
                symmetric = symmetric, call = match.call())
  reslt$svd <- sol
  class(reslt) <- "procrustes"
  reslt
})


