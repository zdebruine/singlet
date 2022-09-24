#' Plot annotations from an NMF model or other compatible objects.
#'
#' @param object      a compatible object (Seurat, DimReduc, nmf, data.frame) 
#' @param ...         additional arguments passed to called functions 
#' @param plot.field  metadata grouping to plot
#' @param reduction   the reduction to plot (default is 'nmf') 
#'
#' @return            a ggplot2 object
#'
#' @export
AnnotationPlot <- function(object, ...) {
  UseMethod("AnnotationPlot")
}


#' @inheritParams AnnotationPlot
#'
#' @rdname AnnotationPlot
#' @name AnnotationPlot
#'
#' @export
#' 
AnnotationPlot.Seurat <- function(object, plot.field = NULL, reduction = "nmf", ...){
  AnnotationPlot(object@reductions[[reduction]])
}


#' @rdname AnnotationPlot
#' @name   AnnotationPlot
#'
#' @export
.S3method("AnnotationPlot", "Seurat", AnnotationPlot.Seurat)


#' Plot metadata enrichment in NMF factors
#' 
#' After running \code{AnnotateNMF}, this function returns 
#' a dot plot of the results
#' 
#' @rdname AnnotationPlot
#' @name AnnotationPlot
#'
#' @inheritParams AnnotationPlot
#'
#' @examples 
#'
#' if (!exists("pbmc3k")) get_pbmc3k_data() %>% NormalizeData -> pbmc3k
#' if (!"nmf" %in% Reductions(pbmc3k)) pbmc3k %>% RunNMF() -> pbmc3k
#' AnnotateNMF(pbmc3k) %>% AnnotationPlot("cell_type")
#'
#' @importFrom stats reshape
#' @importFrom reshape2 melt
#' @import     ggplot2
#'
#' @export
#' 
AnnotationPlot.DimReduc <- function(object, plot.field = NULL, ...){

  if(!("annotations" %in% names(object@misc))){
    stop("the ", reduction, " reduction of this object has no 'annotations' slot. Run 'AnnotateNMF' first.")
  }

  annot <- object@misc$annotations
  if (is.null(plot.field)) {
    plot.field <- names(annot)[[1]]
  } else {
    if(!(plot.field %in% names(annot))) {
      stop(plot.field, "not found in the annotation columns")
    }
    if(length(plot.field) > 1) {
      plot.field <- plot.field[[1]]
    }
  }

  # plot the lods and p-values per factor by group
  AnnotationPlot.data.frame(annot[[plot.field]], plot.field=plot.field)

}


#' @rdname AnnotationPlot
#' @name   AnnotationPlot
#'
#' @export
#'
.S3method("AnnotationPlot", "DimReduc", AnnotationPlot.DimReduc)


#' Plot metadata enrichment in NMF factors
#' 
#' After running \code{AnnotateNMF}, this function returns 
#' a dot plot of the results.  Right now the code is the same as for DimReduc.
#' 
#' @inheritParams AnnotationPlot
#' @rdname AnnotationPlot
#' @name AnnotationPlot
#'
#' @importFrom stats reshape
#' @importFrom reshape2 melt
#' @import     ggplot2
#' @import     RcppML 
#'
#' @export
#' 
AnnotationPlot.nmf <- function(object, plot.field = NULL, ...){

  # nmf objects can have a @misc slot too, so...
  if(!("annotations" %in% names(object@misc))){
    stop("the ", reduction, " reduction of this object has no 'annotations' slot. Run 'AnnotateNMF' first.")
  }

  annot <- object@misc$annotations
  # plot the lods and p-values per factor by group
  AnnotationPlot.data.frame(annot, plot.field=plot.field)

}


#' @rdname AnnotationPlot
#' @name   AnnotationPlot
#'
#' @export
.S3method("AnnotationPlot", "nmf", AnnotationPlot.nmf)


#' Plot metadata enrichment in NMF factors, once summarized into a data.frame
#'
#' @inheritParams AnnotationPlot
#' @rdname AnnotationPlot
#' @name AnnotationPlot
#'
#' @examples 
#' dat <- pbmc3k@reductions$nmf@misc$annotations$cell_type
#' AnnotationPlot(dat, "cell_type")
#'
#' @importFrom stats reshape
#' @importFrom reshape2 melt
#' @import     ggplot2
#'
#' @export
AnnotationPlot.data.frame <- function(object, plot.field, ...) {

  stopifnot(all(c("group", "factor", "fc", "p") %in% names(object)))
  pvals <- reshape(object, timevar = "group", idvar = "factor", 
                   direction = "wide", drop = "fc")
  fc <- reshape(object, timevar = "group", idvar = "factor", 
                direction = "wide", drop = "p") # already adjusted

  rownames(pvals) <- pvals[,1]
  rownames(fc) <- fc[,1]
  pvals <- pvals[, -1]
  fc <- fc[, -1]
  pvals <- as.matrix(pvals)
  fc <- as.matrix(fc)
  colnames(fc) <- colnames(pvals) <- sapply(colnames(pvals), 
                                            function(x) substr(x, 3, nchar(x)))
  fc[fc < 0] <- 0
  idx1 <- hclust(dist(t(fc), method = "binary"), method = "ward.D2")$order
  idx2 <- hclust(dist(fc, method = "binary"), method = "ward.D2")$order

  fc <- fc[idx2, idx1]
  pvals <- pvals[idx2, idx1]
  fc[fc == 0] <- NA
  pvals <- (-1 * log10(pvals))
  pvals[is.infinite(pvals)] <- 100
  pvals[pvals > 100] <- 100

  # these already be melted though?!
  df <- cbind(reshape2::melt(fc), as.vector(pvals))
  colnames(df) <- c("factor", "field", "fc", "pval")
  df$factor <- factor(df$factor, levels = unique(df$factor))

  # construct the plot; return it so the user can tweak it more
  p <- ggplot(df, aes(factor, field, color = pval, size = fc)) + 
         geom_point() + 
         scale_color_viridis_c(option = "B", end = 0.9) + 
         theme_minimal() + 
         labs(y = plot.field, 
              x = "NMF factor", 
              color = "FDR\n(-log10)", 
              size = "association\n(log-odds)") + 
         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
         NULL
  return(p) 

}


#' @rdname AnnotationPlot
#' @name   AnnotationPlot
#'
#' @export
.S3method("AnnotationPlot", "data.frame", AnnotationPlot.data.frame)
