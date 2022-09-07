# singlet v.0.0.99

See the [pkgdown website](https://zdebruine.github.io/singlet/)!

Singlet brings fast Non-negative Matrix Factorization (NMF) with automatic rank determination to the Seurat package for single-cell analysis.

## Install

```{R}
devtools::install_github("zdebruine/singlet")
```

## Vignette

[Guided clustering tutorial with `SeuratData::pbmc3k` dataset](https://zdebruine.github.io/singlet/articles/Guided_Clustering_with_NMF.html)

## Dimension Reduction with NMF

Analyze your single-cell assay with NMF:

```{R}
library(singlet)
library(Seurat)
library(dplyr)
library(cowplot)
pbmc3k <- singlet::get_pbmc3k_data()
set.seed(123)
pbmc3k <- NormalizeData(pbmc3k) %>% 
     RunNMF(k = seq(2, 20, 2), n_replicates = 2) %>%
     RunUMAP(pbmc3k, reduction = "nmf", dims = 1:ncol(pbmc3k@reductions$nmf))

plot_grid(
     RankPlot(pbmc3k) + NoLegend(), 
     DimPlot(pbmc3k) + NoLegend(), 
     labels = "auto", ncol = 2)
```

![NMF workflow](https://github.com/zdebruine/singlet/blob/main/readme_figures/Picture1.png)

NMF can do almost anything that PCA can do, but also:
* imputes missing signal
* always has an optimal rank (for variance-stabilized data)
* uses all the information in your assay (incl. "non-variable" genes)
* is robust across experiments
* learns signatures of transcriptional activity
* is colinear and non-negative (interpretable), rather than orthogonal and signed (not interpretable)

Singlet internally provides the **fastest implementation of NMF**. Cross-validation can take a few minutes for datasets with a few ten thousand cells, but is extremely scalable and runs excellently on HPC nodes and average laptops alike.
