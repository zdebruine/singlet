#' determine the appropriate rank for an ARD ("RADAR") NMF decomposition
#' 
#' @param df  a data.frame of output from crossvalidation: rep, rank, error
#' 
#' @return    the lowest rank that minimizes the reconstruction error
#'
#' @export
GetBestRank <- function(df) {
  if (ncol(df) == 5) {
    # condense to simple format by taking the last iteration in each model
    df <- as.data.frame(group_by(df, rep, k) %>% slice(which.max(iter)))
  }
  if (nrow(df) == 0) {
    return(2)
  } else if (nrow(df) == 1) {
    return(df$k[[1]])
  } else if (length(unique(df$rep)) == 1) {
    return(df$k[[which.min(df$test_error)]])
  } # get the lowest rank for each replicate, take the mean and floor it
  else {
    return(floor(mean((group_by(df, rep) %>% slice(which.min(test_error)))$k)))
  }
}
