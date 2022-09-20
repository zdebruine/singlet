#' @rdname RankPlot
#'
#' @param reduction the NMF reduction slot name (result of \code{RunNMF} where \code{k} was an array)
#' @param detail.level of detail to plot, \code{1} for test set reconstruction error at convergence of each factorization, \code{2} for test set reconstruction error at each fitting iteration of each factorization
#'
#' @export
#'
RankPlot.Seurat <- function(object, reduction = "nmf", detail.level = 1, ...) {
  if (detail.level == 2) {
    plot(subset(object@reductions[[reduction]]@misc$cv_data, iter >= 5), detail.level)
  } else {
    plot(object@reductions[[reduction]]@misc$cv_data, detail.level)
  }
}


#' Plot NMF cross-validation results given a Seurat object
#'
#' S3 method for Seurat that runs the \code{singlet::RunNMF} function.
#'
#' @method RankPlot Seurat
#' @rdname RankPlot
#' @name RankPlot
#'
#' @export
#'
.S3method("RankPlot", "Seurat", RankPlot.Seurat)


#' Plot NMF cross-validation results
#'
#' Given a NMF reduction at multiple ranks, plot rank vs. test set reconstruction error to determine the optimal rank.
#'
#' @param object a Seurat object or a \code{data.frame} that is the result of \code{RunNMF}
#' @param reduction name of the NMF reduction in the Seurat object (result of \code{RunNMF}) for which multiple \code{ranks} were computed.
#' @param ... not implemented
#'
#' @return A ggplot2 object
#'
#' @aliases RankPlot
#'
#' @export
#'
RankPlot <- function(object, reduction = "nmf", ...) {
  UseMethod("RankPlot")
}
