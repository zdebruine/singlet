#' @rdname MetadataSummary
#'
#' @param x result of \code{MetadataSummary}
#'
#' @importFrom reshape2 melt
#'
#' @export
#'
MetadataHeatmap <- function(x) {
  m <- reshape2::melt(as.matrix(x))
  colnames(m) <- c("factor", "group", "frac")
  ggplot(m, aes(x = factor(factor, levels = unique(factor)), y = group, fill = frac)) +
    geom_tile() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank()) +
    labs(x = "factor", y = "group", fill = "relative\ntotal weight") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_gradient2(low = "white", high = "red")
}
