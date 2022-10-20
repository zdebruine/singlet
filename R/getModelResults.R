#' extract data.frame of lods and pvalues for differential factor representation
#' 
#' log-odds of non-null differences for a response by a factor are in fit$lods
#' (which will usually be a matrix), and one-sided p-values for the moderated t
#' test are computed from fit$t and fit$df.total using pt(t, df, lower=FALSE),
#' then adjusted using the step-up procedure of Benjamini & Hochberg.
#' 
#' @param fit   an lmFit result from limma, shrunken with eBayes()
#' @param noneg drop results with negative lods scores? (TRUE) 
#' @param noint drop any results for '(Intercept)'? (TRUE) 
#'
#' @return      a data.frame with columns 'factor', 'group', 'fc', and 'p' 
#' 
#' @details     If an (Intercept) term is found, it will be dropped, and if
#'              negative LODS scores are encountered, they will be dropped,
#'              unless `noneg` and/or `noint` are FALSE.  
#' 
#' @importFrom  reshape2 melt
#' @import      limma
#'
#' @export
getModelResults <- function(fit, noneg=TRUE, noint=TRUE) { 

  # fits are centered, so use signed lods for evidence 
  fcl <- with(fit, melt(lods))
  names(fcl)[3] <- "lods"
  fct <- with(fit, melt(t))
  names(fct)[3] <- "t"
  fcp <- merge(fcl, fct)
  names(fcp)[1:2] <- c("factor", "group")
  fcp$df <- fit$df.total[fcp$factor]
  fcp$p_raw <- with(fcp, pt(t, df, lower=FALSE))
  fcp$p <- NA_real_
    
  # to avoid correcting tests we don't care about,
  # drop (Intercept), negative lods, and p-values
  # (unless requested otherwise) before adjusting
  if (noint) fcp <- subset(fcp, group != "(Intercept)")
  if (noneg) fcp <- subset(fcp, sign(lods) > 0)
  if (nrow(fcp) == 0) {
    message("No associations after filtering.")
  } else { 
    # compute BH-adjusted p-values *after* filtration
    fcp$p <- p.adjust(fcp$p_raw, method="fdr")
  }
  names(fcp) <- sub("^lods$", "fc", names(fcp))
  return(fcp[, c("group", "factor", "fc", "p")])

}
