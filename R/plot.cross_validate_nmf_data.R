#' plot the result of cross-validating rank selection in NMF
#'
#' @param x the result of \code{cross_validate_nmf} (a data.frame)
#' @param detail level of detail to plot
#'
#' @rdname cross_validate_nmf
#'
#' @import ggplot2 
#'
#' @export
#'
plot.cross_validate_nmf_data <- function(x, detail = 2, ...) {
  if (ncol(x) == 5 & detail == 1) {
    x <- as.data.frame(group_by(x, rep, k) %>% slice(which.max(iter)))
    x$iter <- NULL
  }
  if (ncol(x) < 5) {
    x$rep <- factor(x$rep)
    # simple format (detail_level = 1)
    # normalize each replicate to the same minimum
    for (rep in levels(x$rep)) {
      idx <- which(x$rep == rep)
      x$test_error[idx] <- x$test_error[idx] / min(x$test_error[idx])
    }
    best_rank <- GetBestRank(x)
    ggplot(x, aes(k, test_error, color = factor(rep))) +
      geom_point() +
      geom_line() +
      theme_classic() +
      labs(x = "factorization rank", y = "relative test set error", color = "replicate", caption = paste0("(best rank is k = ", best_rank, ")")) +
      theme(aspect.ratio = 1, plot.caption = element_text(hjust = 0.5)) +
      geom_vline(xintercept = best_rank, linetype = "dashed", color = "red") +
      scale_y_continuous(trans = "log10")
  } else {
    # detail_level = 2 format
    best_rank <- GetBestRank(x)
    if (length(unique(x$rep)) == 1) {
      ggplot(x, aes(k, test_error, color = iter, group = iter)) +
        geom_line() +
        scale_color_viridis_c(option = "B") +
        theme_classic() +
        theme(aspect.ratio = 1, plot.caption = element_text(hjust = 0.5)) +
        geom_vline(xintercept = best_rank, linetype = "dashed", color = "red") +
        scale_y_continuous(trans = "log10") +
        labs(x = "factorization rank", y = "test set error", color = "model iteration", caption = paste0("(best rank is k = ", best_rank, ")"))
    } else {
      ggplot(x, aes(k, test_error, color = iter, group = iter)) +
        geom_line() +
        scale_color_viridis_c(option = "B") +
        theme_classic() +
        theme(aspect.ratio = 1, plot.caption = element_text(hjust = 0.5)) +
        geom_vline(xintercept = best_rank, linetype = "dashed", color = "red") +
        scale_y_continuous(trans = "log10") +
        labs(x = "factorization rank", y = "test set error", color = "model iteration", caption = paste0("(best rank is k = ", best_rank, ")")) +
        facet_grid(cols = vars(rep))
    }
  }
}
