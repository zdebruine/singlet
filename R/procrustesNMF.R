#' @title Rescaling of cell and gene features
#'
#' @param Object_X is the reference Object 
#' @param Object_Y is the query object, to be modified to have similiar coordinates.
#' @param scale is used to scale the axes of UMAP coordinates of the Object_Y.
#' @param symmetric to confirm the symmetric of the rotation.
#' @param reduction_X reduction name for the umap reference object.
#' @param reduction_Y query name for the umap reference object.
#' @param reduction.name reduction name of the query object after the procrustes process.
#' @export


setGeneric("procrustesNMF", function(Object_X,Object_Y,...) {
  standardGeneric("procrustesNMF")
})

setMethod("procrustesNMF", signature("Seurat"), function (object_X, object_Y, scale = TRUE, reduction_X="umap_nmf",reduction_Y="umap_nmf_rescaled",
                                                          symmetric = FALSE,reduction.name="umap_nmf_procrustes",...) { 
  X_model <- object_X@reductions[[reduction_X]]
  Y_model <- object_Y@reductions[[reduction_Y]]
  X <- X_model@cell.embeddings
  Y <- Y_model@cell.embeddings
  #X <- scores(X, display = scores, ...)
  #Y <- scores(Y, display = scores, ...)
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
  
  Y_model@cell.embeddings = as.matrix(Yrot)
  object_Y@reductions[[reduction.name]] <- Y_model
  object_Y
})

