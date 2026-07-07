#' Plot GSEA results on a heatmap
#'
#' Plot top GSEA terms for each NMF factor on a heatmap.
#'
#' @param object A Seurat, SingleCellExperiment, or RcppML::nmf object,
#'   a list of GSEA results, or a matrix of -log10 adjusted p-values.
#' @param reduction Dimensional reduction for which GSEA analysis has been performed.
#' @param collapse collapse slot in GSEA result
#' @param padj.sig Adjusted p-value threshold for selecting significant pathways.
#' @param max.terms.per.factor Maximum number of top pathways to show per NMF factor.
#' @param dropcommon Whether to remove pathways significant in more than half
#'   of the NMF factors.
#' @param gsea.name Name under which GSEA results are stored. Default is "gsea".
#' @param truncate.terms Maximum number of characters to retain in pathway names.
#'   Set to NULL to keep full names.
#' @param ... Additional arguments.
#'
#' @return A ggplot2 object.
#'
#' @export
GSEAHeatmap <- function(object, ...) {
  UseMethod("GSEAHeatmap")
}

#' @rdname GSEAHeatmap
#' @method GSEAHeatmap Seurat
#' @export
GSEAHeatmap.Seurat <- function(object, reduction = "nmf",
                               collapse = NULL,
                               padj.sig = 0.05,
                               max.terms.per.factor = 3,
                               dropcommon = TRUE,
                               gsea.name = "gsea",
                               truncate.terms = 48,
                               ...) {
  gsea <- object@reductions[[reduction]]@misc[[gsea.name]]
  
  if (is.null(gsea)) {
    stop("GSEA results not found. Run RunGSEA() first.")
  }
  
  p <- GSEAHeatmap(
    object = gsea,
    collapse = collapse,
    padj.sig = padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  
  return(p)
}
.S3method("GSEAHeatmap", "Seurat", GSEAHeatmap.Seurat)

#' @rdname GSEAHeatmap
#' @method GSEAHeatmap SingleCellExperiment
#' @export
GSEAHeatmap.SingleCellExperiment <- function(object, reduction = "nmf",
                                             collapse = NULL,
                                             padj.sig = 0.05,
                                             max.terms.per.factor = 3,
                                             dropcommon = TRUE,
                                             gsea.name = "gsea",
                                             truncate.terms = 48,
                                             ...) {
  reduc = object@metadata[[reduction]]
  
  p <- GSEAHeatmap(
    object = reduc,
    collapse = collapse,
    padj.sig = padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    gsea.name = gsea.name,
    truncate.terms = truncate.terms
  )
  
  return(p)
}

.S3method("GSEAHeatmap", "SingleCellExperiment", GSEAHeatmap.SingleCellExperiment)

#' @rdname GSEAHeatmap
#' @method GSEAHeatmap nmfList
#' @export
GSEAHeatmap.nmfList <- function(object,
                                collapse = NULL,
                                padj.sig = 0.05,
                                max.terms.per.factor = 3,
                                dropcommon = TRUE,
                                gsea.name = "gsea",
                                truncate.terms = 48,
                                ...) {
  
  if (is.null(object[[gsea.name]])) {
    stop("GSEA results not found. Run RunGSEA() first.")
  }
  p <- GSEAHeatmap(
    object = object[[gsea.name]],
    collapse = collapse,
    padj.sig = padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  return(p)
}

.S3method("GSEAHeatmap", "nmfList", GSEAHeatmap.nmfList)

#' @rdname GSEAHeatmap
#' @method GSEAHeatmap nmf
#' @export
GSEAHeatmap.nmf <- function(object,
                            collapse = NULL,
                            padj.sig = 0.05,
                            max.terms.per.factor = 3,
                            dropcommon = TRUE,
                            gsea.name = "gsea",
                            truncate.terms = 48,
                            ...) {
  
  if (is.null(object@misc[[gsea.name]])) {
    stop("GSEA results not found. Run RunGSEA() first.")
  }
  p <- GSEAHeatmap(
    object = object@misc[[gsea.name]],
    collapse = collapse,
    padj.sig = padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  
  return(p)
}

.S3method("GSEAHeatmap", "nmf", GSEAHeatmap.nmf)

#' @rdname GSEAHeatmap
#' @method GSEAHeatmap list
#' @importFrom fgsea collapsePathways
#' @export
GSEAHeatmap.list <- function(object,
                             collapse = NULL,
                             padj.sig = 0.05,
                             max.terms.per.factor = 3,
                             dropcommon = TRUE,
                             truncate.terms = 48,
                             ...) {
  if (!all(c("gsea", "ranks", "pathways") %in% names(object))) {
    stop("GSEA results not found. Run RunGSEA() first.")
  }
  
  if(!is.null(collapse)){
    if(!(collapse %in% names(object))){
      stop("Collapse slot does not exist!")
      }
    collapsed =  object[[collapse]]
    main_pathways = getMainPathways(collapsed)
    formatted =  formatGSEA(
      object,
      pathways = main_pathways
    )
  } else{
    formatted = formatGSEA(object)
  }
  
  p <- GSEAHeatmap(
    object = formatted[["padj"]],
    padj.sig = padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  
  return(p)
}

.S3method("GSEAHeatmap", "list", GSEAHeatmap.list)



#' @rdname GSEAHeatmap
#' @method GSEAHeatmap matrix
#' @export
GSEAHeatmap.matrix <- function(object,
                               padj.sig = 0.05,
                               max.terms.per.factor = 3,
                               dropcommon = TRUE,
                               truncate.terms = 48,
                               ...) {
  df <- selectGSEATerms(
    df = object,
    padj.sig=padj.sig,
    max.terms.per.factor = max.terms.per.factor,
    dropcommon = dropcommon,
    truncate.terms = truncate.terms
  )
  
  if (nrow(df) == 0) {
    stop("No significant GSEA pathways remain after filtering")
  }
  
  df <- reshape2::melt(df)
  p <- plotGSEAHeatmap(df)
  
  return(p)
}

.S3method("GSEAHeatmap", "matrix", GSEAHeatmap.matrix)


plotGSEAHeatmap = function(df){
  p <- ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "B") +
    theme_classic() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      x = "NMF factor",
      y = "Gene Set",
      fill = "FDR\n(-log10)"
    ) +
    theme(
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) + 
    NULL 
  
  return(p) 
}

selectGSEATerms = function(df, padj.sig = 0.05,max.terms.per.factor = 3,dropcommon = TRUE,truncate.terms=48){
  df = df[(rowSums(df,na.rm=TRUE) > 0 ), , drop = FALSE]
  if (nrow(df) == 0) return(df)
  # markers for each factor based on the proportion of signal in that factor
  df2 <- as.matrix(Diagonal(x = 1 / rowSums(df,na.rm = TRUE)) %*% df)
  
  # see https://github.com/zdebruine/singlet/issues/26
  # thanks to @earbebarnes
  rownames(df2) <- rownames(df) #add row names to df2
  
  terms <- c()
  for (i in seq_len(ncol(df2))) {
    terms_i <- df[, i]
    idx <- terms_i > -log10(padj.sig)
    terms_i <- terms_i[idx]
    terms_j <- df2[idx, i]
    v <- sort(terms_j, decreasing = TRUE)
    v <- v[!is.na(v)] #making this explicit
    if (length(v) > max.terms.per.factor) {
      terms <- c(terms, names(v)[1:max.terms.per.factor])
    } else {
      terms <- c(terms, names(v))
    }
  }
  terms <- unique(terms)
  df <- df[terms, , drop = FALSE]
  #could be earlier, but this lets us exclude NAs where 0 could be < max.terms 
  df[is.na(df)] = 0
  
  if(!is.null(truncate.terms)){
    rownames(df) <- sapply(rownames(df), function(x) {
      ifelse(nchar(x) > truncate.terms, paste0(substr(x, 1, truncate.terms), "..."), x)
    })
  }
  
  if (dropcommon && ncol(df) > 1) { 
    # remove terms that are broadly significant
    v <- which((rowSums(df > -log10(padj.sig)) > (ncol(df) / 2)))
    if (length(v) > 0) df <- df[-v, , drop = FALSE]
  }
  
  return(df)
}