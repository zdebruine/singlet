#' @inheritParams AnnotateNMF.Seurat
#' @aliases AnnotationPlot
#' @rdname AnnotateNMF
#' @export
#' 
AnnotationPlot.Seurat <- function(object, plot.field = NULL, reduction = "nmf", ...){
  AnnotationPlot(object@reductions[[reduction]])
}


#' @export
#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
.S3method("AnnotationPlot", "Seurat", AnnotationPlot.Seurat)


#' Plot metadata enrichment in NMF factors
#' 
#' After running \code{AnnotateNMF}, this function returns 
#' a dot plot of the results
#' 
#' @inheritParams AnnotateNMF.DimReduc
#' @param plot.field name of field in \code{meta.data} to plot, if \code{NULL}, plots the first field in the annotation results
#' @return ggplot2 object
#' @aliases AnnotationPlot
#' @rdname AnnotateNMF
#' @export
#' @importFrom stats reshape t.test
#' 
AnnotationPlot.DimReduc <- function(object, plot.field = NULL, ...){

  if(!("annotations" %in% names(object@misc))){
    stop("the ", reduction, " reduction of this object has no 'annotations' slot. Run 'AnnotateNMF' first.")
  }

  if(is.null(field)){
    field <- names(object@misc$annotations)[[1]]
  } else {
    if(!(field %in% names(object@misc$annotations))){
      stop("specified field was not in the annotation object")
    }
    if(length(field) > 1) field <- field[[1]]
  }

  # construct a matrix of pvalues and fc
  ann <- object@misc$annotations[[field]]

  pvals <- reshape(ann, timevar = "group", 
                   idvar = "factor", direction = "wide", drop = "fc")
  fc <- reshape(ann, timevar = "group", 
                idvar = "factor", direction = "wide", drop = "p")
  rownames(pvals) <- pvals[,1]
  rownames(fc) <- fc[,1]
  pvals <- pvals[, -1]
  fc <- fc[, -1]
  pvals <- as.matrix(pvals)
  fc <- as.matrix(fc)
  colnames(fc) <- colnames(pvals) <- sapply(colnames(pvals), 
                                            function(x) substr(x, 3, nchar(x)))
  fc[fc < 0] <- 0
  idx1 <- hclust(dist(t(fc)), method = "ward.D2")$order
  idx2 <- hclust(dist(fc), method = "ward.D2")$order
  fc <- fc[idx2, idx1]
  pvals <- pvals[idx2, idx1]
  fc[fc == 0] <- NA
  pvals <- -log10(pvals)
  pvals[is.infinite(pvals)] <- 100
  pvals[pvals > 100] <- 100
  df <- cbind(reshape2::melt(fc), as.vector(pvals))
  colnames(df) <- c("factor", "field", "fc", "pval")
  df$factor <- factor(df$factor, levels = unique(df$factor))
  p <- ggplot(df, 
              aes(factor, field, color = pval, size = fc)) + 
         geom_point() + 
         scale_color_viridis_c(option = "B", end = 0.9) + 
         theme_classic() + 
         labs(y = field, 
              x = "NMF factor", 
              color = "p-value\n(-log10)", 
              size = "fold enrichment") + 
         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
         NULL

  # can just return p and print it as needed
  # users will have it autoplotted,
  # and can modify the plot object 
  # suppressWarnings(print(p))
  return(p) 

}


#' @export
#' @rdname AnnotateNMF
#' @name AnnotateNMF
#'
.S3method("AnnotationPlot", "DimReduc", AnnotationPlot.DimReduc)


#' @export
#' @rdname AnnotateNMF
#'
AnnotationPlot <- function(object, ...) {
  UseMethod("AnnotationPlot")
}


