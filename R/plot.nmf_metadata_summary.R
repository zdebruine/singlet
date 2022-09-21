#' @rdname MetadataSummary
#'
#' @param x   a data.frame
#' @param ... not implemented
#'
#' @importFrom reshape2 melt
#'
#' @export
#'
plot.nmf_metadata_summary <- function(x, ...) {
  m <- reshape2::melt(as.matrix(x))
  colnames(m) <- c("group", "factor", "frac")
  ggplot(m, aes(x = factor(factor, levels = unique(factor)), y = frac, fill = group)) +
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "factor", y = "Representation in group") +
    scale_y_continuous(expand = c(0, 0))
}


#' @rdname MetadataSummary
#'
#' @name MetadataSummary
#'
#' @export
#'
.S3method("plot", "nmf_metadata_summary", plot.nmf_metadata_summary)
