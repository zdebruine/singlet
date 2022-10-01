# Row-wise rasterization of a sparse matrix
#' 
#' Bin together values from every block of \code{n} rows and calculate mean value, with a sparse \code{dgCMatrix} as input and a dense \code{matrix} as output. This technique is useful in some genomics applications.
#'
#' @param A matrix to be rasterized
#' @param n row-wise binning size
#' @param threads number of threads to use (0 to let OpenMP decide how many are available and use them all)
#' @export
#' 
RasterizeRowwise <- function(A, n = 10, threads = 0){
  B <- rowwise_compress(A, n, threads)
  rownames(B) <- rownames(A)[seq(1, floor(nrow(A) / n) * n, n)]
  colnames(B) <- colnames(A)
  B
}