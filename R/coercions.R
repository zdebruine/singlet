#' @exportMethod coerce
#' @importClassesFrom RcppML nmf
setAs("list", "nmf", 
      function(from) {

        msg <- NULL
        required <- c("w", "d", "h")
        if (!all(required %in% names(from))) {
          msg <- c(msg, 
            "Cannot find $w, $d, and $h to create nmf object from list.")
        } else { 
          if (ncol(from$w) != nrow(from$h)) {
            msg <- c(msg, 
              "The $w and $h matrices are of unequal rank. Cannot coerce.")
          }
          if (ncol(from$w) != length(from$d)) {
            msg <- c(msg, 
              "The scaling diagonal $d is the wrong length. Cannot coerce.")
          }
        }
      
        if (!is.null(msg)) {
          stop(msg)
        } else { 
          new("nmf",
              w = from$w,
              d = from$d,
              h = from$h,
              misc = from[setdiff(names(from), required)])
        }

      })


#' @exportMethod coerce
#' @importClassesFrom RcppML nmf
if (requireNamespace("SingleCellExperiment", quietly=TRUE)) {
  setAs("nmf", "LinearEmbeddingMatrix", function(from) {
    factorNames <- colnames(from@w)
    sampleNames <- colnames(from@h) 
    lem <- LinearEmbeddingMatrix(sampleFactors=t(from@h), 
                                 featureLoadings=from@w,
                                 factorData=DataFrame(d=from@d, 
                                                      row.names=factorNames),
                                 metadata=from@misc)
    rownames(lem) <- sampleNames
    return(lem)
  })
}


#' @exportMethod coerce
#' @importClassesFrom RcppML nmf
if (requireNamespace("SingleCellExperiment", quietly=TRUE)) {
  setAs("LinearEmbeddingMatrix", "nmf", function(from) {
    d <- factorData(from)$d
    names(d) <- rownames(factorData(from))
    new("nmf", 
        w = featureLoadings(from), 
        d = d,
        h = t(sampleFactors(from)),
        misc = metadata(from))
  })
}
