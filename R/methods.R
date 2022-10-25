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


#' @exportMethod [ 
#' @importClassesFrom RcppML nmf
setMethod("[", "nmf",
          function(x, i, j, ..., drop=TRUE) {
            if (missing(i) & missing(j)) return(x)
            if (missing(i)) i <- seq_along(x@d)
            if (missing(j)) j <- colnames(x@h)
            if ("covs" %in% names(x@misc)) x@misc$covs <- x@misc$covs[j, ]
            new("nmf", w = x@w[, i], d = x@d[i], h = x@h[i, j], misc = x@misc)
          })


#' @exportMethod $
#' @importClassesFrom RcppML nmf
setMethod("$", "nmf", 
          function(x, name) {
            if ("covs" %in% names(x@misc)) {
              x@misc$covs[[name]]
            } else { 
              NULL
            }
          })


#' @exportMethod $<-
#' @importClassesFrom RcppML nmf
setReplaceMethod("$", "nmf", 
          function(x, name, value) {
            if (is.null(x@misc$covs)) {
              x@misc$covs <- data.frame(row.names = colnames(x@h))
            }
            x@misc$covs[[name]] <- value
            return(x)
          })
