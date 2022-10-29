#' determine the appropriate rank for an AutoNMF decomposition
#'
#' @param df  a data.frame of output from crossvalidation: rep, rank, error
#' @inheritParams RunNMF
#' @return    the lowest rank that minimizes the reconstruction error
#' @export
#'
GetBestRank <- function(df, tol.overfit = 1e-4, ...) {
  df$rep <- factor(df$rep)
  best_ranks <- c()
  for (replicate in levels(df$rep)) {
    df_rep <- subset(df, rep == replicate)
    # calculate overfitting tolerance
    max_rank <- max(df_rep$k) + 1
    for (rank in unique(df_rep$k)) {
      if (rank < max_rank) {
        df_rank <- subset(df_rep, k == rank)
        if (nrow(df_rank) > 1) {
          v2 <- df_rank$test_error[2:nrow(df_rank)]
          v1 <- df_rank$test_error[1:(nrow(df_rank) - 1)]
          for (pos in 2:length(v1)) {
            if (v1[[pos]] > v1[[pos - 1]]) v1[[pos]] <- v1[[pos - 1]]
          }
          if (max(c(0, (v2 - v1) / (v2 + v1))) > tol.overfit) {
            max_rank <- rank
          }
        }
      }
    }
    df_rep <- subset(df_rep, k < max_rank)
    if (nrow(df_rep) == 0) {
      best_ranks <- c(best_ranks, 2)
    } else if (nrow(df) == 1) {
      best_ranks <- c(best_ranks, df_rep$k[[1]])
    } else {
      # condense to simple format by taking the last iteration in each model
      df_rep <- as.data.frame(group_by(df_rep, rep, k) %>% slice(which.max(iter)))
      best_ranks <- c(best_ranks, df_rep$k[which.min(df_rep$test_error)])
    }
  }

  # get the lowest rank for each replicate, take the mean and floor it
  floor(mean(best_ranks))
}
