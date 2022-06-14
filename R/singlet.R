#' Singlet
#' 
#' Fast single-cell analysis with non-negative dimensional reductions
#' 
#' @details 
#' There are reasons to not use PCA. 
#' * PCA fits to missing signal, 
#' * considers only highly variable features, 
#' * is almost useless without further graph-based analysis, 
#' * requires centering and scaling of your data,
#' * and is robust only within experiments.
#' 
#' Instead, you should use Non-negative Matrix Factorization (NMF).
#' * NMF imputes missing signal,
#' * learns models using all features,
#' * does everything PCA does and provides useful information itself,
#' * requires only variance stabilization,
#' * and is robust across experiments.
#' 
#' Singlet is all about extremely fast NMF for single-cell dimensional reduction and integration.
#' 
#' See the vignettes to get started.
#' 
#' @import ggplot2 Seurat dplyr reshape2 Matrix utils stats methods msigdbr fgsea
#' @useDynLib singlet, .registration = TRUE
#' @docType package
#' @name singlet
#' @author Zach DeBruine
#' @aliases singlet-package
#' @md
#'
NULL