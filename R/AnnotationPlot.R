#' Plot annotations from an NMF model or other compatible objects.
#'
#' @param object      a compatible object (Seurat, DimReduc, nmf, data.frame) 
#' @param ...         additional arguments passed to called functions 
#' @param plot.field  metadata grouping to plot
#' @param reduction   the reduction to plot (default is 'nmf') 
#' @param dropEmpty   drop factors without significant associations? (TRUE)
#' @param annotation.name name of factor annotation, default = "annotations"
#'
#' @return            a ggplot2 object
#'
#' @export
AnnotationPlot <- function(object, ...) {
  UseMethod("AnnotationPlot")
}


#' @rdname AnnotationPlot
#' @name AnnotationPlot
#'
#' @export
#' 
AnnotationPlot.Seurat <- function(object, plot.field = NULL, reduction = "nmf", dropEmpty=TRUE, annotation.name = "annotations", cluster_by="pos", ...){
  AnnotationPlot(object@reductions[[reduction]],plot.field = plot.field, dropEmpty=dropEmpty, annotation.name = annotation.name, cluster_by = cluster_by)
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
#' @examples 
#' \dontrun{
#' get_pbmc3k_data() %>% NormalizeData %>% RunNMF %>% AnnotateNMF -> pbmc3k
#' AnnotationPlot(pbmc3k, "cell_type")
#' }
#' @importFrom stats reshape
#' @importFrom reshape2 melt
#' @import     ggplot2
#'
#' @export
#' 
AnnotationPlot.DimReduc <- function(object, plot.field=NULL, dropEmpty=TRUE, annotation.name = "annotations", cluster_by="pos", ...) {
  
  if(!(annotation.name %in% names(object@misc))){
    stop("the ", reduction, " reduction of this object has no 'annotations' slot. Run 'AnnotateNMF' first.")
  }
  
  annot <- object@misc[[annotation.name]]
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
  
  # plot per factor by group
  AnnotationPlot.data.frame(annot[[plot.field]], 
                            plot.field=plot.field,
                            dropEmpty=dropEmpty,
                            cluster_by = cluster_by)
  
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
AnnotationPlot.nmf <- function(object, plot.field=NULL, dropEmpty=TRUE, annotation.name = "annotations", cluster_by="pos",...) {
  
  # nmf objects can have a @misc slot too, so...
  if(!("annotations" %in% names(object@misc))){
    stop("the ", reduction, " reduction of this object has no 'annotations' slot. Run 'AnnotateNMF' first.")
  }
  
  annot <- object@misc[[annotation.name]]
  AnnotationPlot(annot, plot.field=plot.field, dropEmpty=dropEmpty,cluster_by = cluster_by)
  
}


#' @rdname AnnotationPlot
#' @name   AnnotationPlot
#'
#' @export
.S3method("AnnotationPlot", "nmf", AnnotationPlot.nmf)


#' Plot metadata enrichment in NMF factors, if a list of data.frames
#'
#' @rdname AnnotationPlot
#' @name AnnotationPlot
#'
#' @export
AnnotationPlot.list <- function(object, plot.field, dropEmpty=TRUE, cluster_by="pos",...) {
  
  stopifnot(plot.field %in% names(object))
  AnnotationPlot.data.frame(object[[plot.field]], 
                            plot.field=plot.field, 
                            dropEmpty=dropEmpty,
                            cluster_by = cluster_by,
                            ...)
  
}


#' @rdname AnnotationPlot
#' @name   AnnotationPlot
#'
#' @export
.S3method("AnnotationPlot", "list", AnnotationPlot.list)



#' Plot metadata enrichment in NMF factors, once summarized into a data.frame
#'
#' @rdname AnnotationPlot
#' @name AnnotationPlot
#'
#' @examples 
#' \dontrun{
#' dat <- pbmc3k@reductions$nmf@misc$annotations$cell_type
#' AnnotationPlot(dat, "cell_type")
#'
#' # if running interactively:
#' library(plotly)
#' ggplotly(AnnotationPlot(dat, "cell_type"))
#' }
#' @importFrom stats reshape
#' @importFrom reshape2 melt
#' @importFrom reshape2 acast
#' @import     ggplot2
#'
#' @export
AnnotationPlot.data.frame <- function(object, plot.field, dropEmpty=TRUE, cluster_by="pos",...){
  
  pcols <- c("factor","group","p")
  fcols <- c("factor","group","fc")
  ccols = c("factor","group","coefficients")
  mandatory <- union(union(pcols, fcols),ccols)
  stopifnot(all(mandatory %in% names(object)))
  
  # cast and cluster LODS ("fc"); could stand to be sparse
  fc <- reshape2::acast(object[, fcols], group ~ factor, value.var="fc")
  pvals <- reshape2::acast(object[, pcols], group ~ factor, value.var="p")
  coeff <- reshape2::acast(object[, ccols], group ~ factor, value.var="coefficients")
  # drop "nmf" from the names
  # colnames(fc) <- sub("^nmf", "", colnames(fc)) 
  # colnames(pvals) <- sub("^nmf", "", colnames(pvals)) 
  
  # check that all dimnames match
  stopifnot(identical(rownames(fc), rownames(pvals)))
  stopifnot(identical(colnames(fc), colnames(pvals)))
  
  # clean up for clustering on LODs ("fold change")
  fc[fc < 0] <- 0      # drop out factors with no enrichment
  fc[is.na(fc)] <- 0   # replace NAs so that clustering works
  fdr_weight <- round(-1 * log10(pvals)) # cut off at an FDR of ~0.317
  fc[fdr_weight == 0] <- 0 # fdr > 0.317
  
  # cluster on "fold change" (LODS)
  # ridx <- rev(hclust(dist(fc, method = "binary"), method = "ward.D2")$order)
  # fields <- rownames(fc)[ridx]
  # cidx <- rev(hclust(dist(t(fc), method = "binary"), method = "ward.D2")$order)
  # factors <- colnames(fc)[cidx]
  fields = rownames(fc)[getClusterIdx(fc,coeff,cluster_by)]
  factors = colnames(fc)[getClusterIdx(t(fc),t(coeff),cluster_by)]
  
  pvals <- pvals[fields, factors]
  fc <- fc[fields, factors]
  fc[fc == 0] <- NA 
  
  # check that all dimnames match
  stopifnot(identical(rownames(fc), rownames(pvals)))
  stopifnot(identical(colnames(fc), colnames(pvals)))
  
  # rescale adjusted p-values to -1 * log10(fdr)
  negative_log10_fdr <- round(-1 * log10(pvals))
  negative_log10_fdr[is.infinite(negative_log10_fdr)] <- 100
  negative_log10_fdr[negative_log10_fdr > 100] <- 100
  is.na(negative_log10_fdr[which(negative_log10_fdr == 0)]) <- TRUE
  
  # check (again!) that all dimnames match
  stopifnot(identical(rownames(fc), rownames(negative_log10_fdr)))
  stopifnot(identical(colnames(fc), colnames(negative_log10_fdr)))
  
  # melt both and merge
  df <- merge(
    merge(melt(fc, value.name="lods"), 
          melt(negative_log10_fdr, value.name="negative_log10_fdr")),
    melt(coeff, value.name="coefficients")
  )
  names(df) <- c("field", "factor", "lods", "negative_log10_fdr","coefficients")
  
  # retain clustering of NMF fractors along the columns
  df$factor <- factor(df$factor, levels=factors)
  
  # retain clustering of predictors along the rows
  df$field <- factor(df$field, levels=fields)
  
  # drop factors and fields without associations, unless requested not to 
  if (dropEmpty) df <- subset(df, !is.na(lods) & !is.na(negative_log10_fdr))
  
  # not always useful:
  df$design <- plot.field
  
  # condense somewhat 
  df <- df[, c("design", "field", "factor", "lods", "negative_log10_fdr","coefficients")]
  df <- df[rev(order(df$negative_log10_fdr)), ]
  
  #caution: multiplying by the sign of coefficients will only work if coefficients have no 0
  #this may look more confusing if the FDR contain non significant results, which may already be implicitly filtered.
  # construct plot
  p <- ggplot(df, 
              aes(x = factor, 
                  y = field, 
                  color = negative_log10_fdr*sign(coefficients), 
                  fill = negative_log10_fdr*sign(coefficients),
                  size = lods,
                  shape = factor(ifelse(coefficients > 0, "Positive", "Negative"),
                                 levels = c("Positive", "Negative"))
              )
  ) + 
    geom_point() + 
    scale_shape_manual(values = c("Positive" = 24, "Negative" = 25)) +
    scale_colour_gradientn(
      colours = c("blue", "white", "red"),
      values = scales::rescale(c(-100, -50 , log10(0.05) , 0, -log10(0.05) , 50 , 100)),
      limits = c(-100, 100),
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar"
    )+
    scale_fill_gradientn(
      colours = c("blue", "white", "red"),
      values = scales::rescale(c(-100, -50 , log10(0.05), 0, -log10(0.05), 50 ,100)),
      limits = c(-100, 100),
      space = "Lab",
      na.value = "grey50",
      guide = "none"
    )+
    # scale_color_gradient2(
    #   low = "blue", mid = "white", high = "red", midpoint = 0
    # )+
    # scale_fill_gradient2(
    #   low = "blue", mid = "white", high = "red", midpoint = 0
    # ,guide = "none")+
    # scale_color_viridis_c(direction = -1, option = "B", end = 0.9) +
    # scale_fill_viridis_c(direction = -1, option = "B", end = 0.9,guide = "none") +
    theme_minimal() + 
    labs(y = plot.field, 
         x = "Factor",
         shape = "Direction",
         color = " sign(coeffs) x\n-log10(FDR)", 
         size = "Association\n(log-odds)"
    ) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    guides(alpha = "none") + 
    NULL
  
  # return so the user can tweak it if necessary
  return(p) 
  
}


#' @rdname AnnotationPlot
#' @name   AnnotationPlot
#'
#' @export
.S3method("AnnotationPlot", "data.frame", AnnotationPlot.data.frame)

getClusterIdx = function(fc,coeff,cluster_by="pos"){
  stopifnot(cluster_by %in% c("std","pos","neg"))
  
  if(cluster_by == "pos"){
    fc[sign(coeff)==-1] = 0
  }else if (cluster_by=="neg"){
    fc[sign(coeff)==1] = 0
  }
  
  idx = rev(hclust(dist(fc, method = "binary"), method = "ward.D2")$order)
  
  return(idx)
}