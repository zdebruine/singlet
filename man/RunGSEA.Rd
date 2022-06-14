% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_wrapper.R
\name{RunGSEA}
\alias{RunGSEA}
\title{Run Gene Set Enrichment Analysis on a Reduction}
\usage{
RunGSEA(
  object,
  reduction = "nmf",
  species = "Homo sapiens",
  category = "C5",
  min.size = 25,
  max.size = 500,
  collapse = TRUE,
  dims = NULL,
  verbose = TRUE,
  padj.sig = 0.01,
  ...
)
}
\arguments{
\item{object}{a Seurat object}

\item{reduction}{dimensional reduction to use}

\item{species}{species for which to load gene sets}

\item{category}{msigdbr gene set category (i.e. "H", "C5", etc.)}

\item{min.size}{minimum number of terms in a gene set}

\item{max.size}{maximum number of terms in a gene set}

\item{collapse}{filter pathways to remove highly redundant terms}

\item{dims}{factors in the reduction to use, default \code{NULL} for all factors}

\item{verbose}{print progress to console}

\item{...}{additional arguments to \code{fgseaMultilevel}}

\item{padj.size}{significance cutoff for BH-adjusted p-values (default 0.01)}
}
\value{
a Seurat object, with GSEA information in the misc slot. BH-adj p-values are on a -log10 scale.
}
\description{
Run GSEA to identify gene sets that are enriched within NMF factors.
}