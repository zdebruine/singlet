#' Summarize contribution of sample groups to NMF factors
#'
#' Calculate the mean weight of samples in discrete and unique groups to each factor
#'
#' @rdname MetadataSummary
#'
#' @param h matrix giving factors as rows and samples as columns
#' @param factor_data a factor of the same length as the number of columns in \code{h}
#' @param reorder sort results by proportion in each group (uses \code{hclust} if >2 groups)
#'
#' @return \code{data.frame} of mean weights for each sample group within each factor of class \code{nmf_metadata_summary}. Use the \code{plot} method to visualize.
#'
#' @export
#'
MetadataSummary <- function(h, factor_data, reorder = TRUE) {
  factor_data <- as.factor(factor_data)
  if (is.null(rownames(h))) rownames(h) <- paste0("factor", 1:nrow(h))
  m <- matrix(0, nrow(h), length(levels(factor_data)))
  rownames(m) <- rownames(h)
  colnames(m) <- levels(factor_data)
  for (j in 1:length(levels(factor_data))) {
    for (i in 1:nrow(h)) {
      m[i, j] <- mean(h[i, which(factor_data == levels(factor_data)[[j]])])
    }
  }
  m <- apply(m, 1, function(x) x / sum(x))
  if (length(levels(factor_data)) == 2) {
    m <- m[order(m[, 1], decreasing = TRUE), ]
  } else if (reorder) {
    m <- m[hclust(dist(m), method = "ward.D2")$order, hclust(dist(t(m)), method = "ward.D2")$order]
  }
  t(m)
  m <- as.data.frame(m)
  class(m) <- c("nmf_metadata_summary", "data.frame")
  m
}
