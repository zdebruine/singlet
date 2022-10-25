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
