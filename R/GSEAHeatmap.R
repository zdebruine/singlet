#' Plot GSEA results on a heatmap
#'
#' Plot top GSEA terms for each NMF factor on a heatmap
#'
#' @param object Seurat or RcppML::nmf object
#' @param reduction a dimensional reduction for which GSEA analysis has been performed
#' @param max.terms.per.factor show this number of top terms for each factor
#' @param dropcommon  drop broadly enriched terms across factors? (TRUE) 
#'
#' @return ggplot2 object
#'
#' @export
#'
GSEAHeatmap <- function(object, reduction = "nmf", max.terms.per.factor = 3, dropcommon = TRUE) {

  if (is(object, "Seurat")) {
    df <- object@reductions[[reduction]]@misc$gsea$padj
  } else if (is(object, "nmf")) {
    df <- object@misc$gsea$padj
  }
  
  # markers for each factor based on the proportion of signal in that factor
  df2 <- as.matrix(Diagonal(x = 1 / rowSums(df)) %*% df)

  # see https://github.com/zdebruine/singlet/issues/26
  # thanks to @earbebarnes
  rownames(df2) <- rownames(df) #add row names to df2

  terms <- c()
  for (i in 1:ncol(df2)) {
    terms_i <- df[, i]
    idx <- terms_i > -log10(0.05)
    terms_i <- terms_i[idx]
    terms_j <- df2[idx, i]
    v <- sort(terms_j, decreasing = TRUE)
    if (length(v) > max.terms.per.factor) {
      terms <- c(terms, names(v)[1:max.terms.per.factor])
    } else {
      terms <- c(terms, names(v))
    }
  }
  terms <- unique(terms)
  df <- df[terms, ]

  rownames(df) <- sapply(rownames(df), function(x) {
    ifelse(nchar(x) > 48, paste0(substr(x, 1, 45), "..."), x)
  })

  if (dropcommon) { 
    # remove terms that are broadly significant
    v <- which((rowSums(df > -log10(0.05)) > (ncol(df) / 2)))
    if (length(v) > 0) df <- df[-v, ]
  }
  df <- reshape2::melt(df)
  p <- ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "B") +
    theme_classic() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      x = "NMF factor",
      y = "GO Term",
      fill = "FDR\n(-log10)"
    ) +
    theme(
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) + 
    NULL 

  return(p) 

}
