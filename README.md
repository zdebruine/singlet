# singlet v0.99.9

See the [pkgdown website](https://zdebruine.github.io/singlet/)!

`singlet` brings fast Non-negative Matrix Factorization (NMF) with automatic rank determination to Seurat and
SingleCellExperiment objects.

**Status.** This package is being folded into the [singlet.bio](https://singlet.bio) project
([Singlet-Bio/singlet](https://github.com/Singlet-Bio/singlet)), which loads the atlas of reprocessed public
single-cell studies and will carry these functions forward under the same names. Until then this version is
maintained: it works with Seurat v3/v4 (`Assay`) and v5 (`Assay5`) objects, and the rank search no longer crashes on
rank 1 or on matrices with no signal. See `docs/` in the Singlet-Bio repository for the migration plan.

## Install

Install the Bioconductor dependencies, then `singlet`:

```{R}
install.packages(c("BiocManager", "devtools"))
BiocManager::install(c("fgsea", "limma", "scuttle", "SingleCellExperiment"))
devtools::install_github("zdebruine/singlet")
```

`singlet` builds against `RcppML` (CRAN, version 1.0.0 was used for the tests) and needs a C++17 compiler: Rtools on
Windows, Xcode command line tools on macOS, `build-essential` on Linux.

## Introductory Vignette

[Guided clustering tutorial](https://zdebruine.github.io/singlet/articles/Guided_Clustering_with_NMF.html)

## Dimension Reduction with NMF

Analyze your single-cell assay with NMF:

```{R}
library(singlet)
library(Seurat)
library(dplyr)
library(cowplot)
set.seed(123) # for reproducible NMF models
get_pbmc3k_data() %>% NormalizeData %>% RunNMF -> pbmc3k
pbmc3k <- RunUMAP(pbmc3k, reduction = "nmf", dims = 1:ncol(pbmc3k@reductions$nmf))

plot_grid(
     RankPlot(pbmc3k) + NoLegend(), 
     DimPlot(pbmc3k) + NoLegend(), 
     ncol = 2)
```

NMF can do almost anything that PCA can do, but also imputes missing signal, always has an optimal rank (for variance-stabilized data), uses all the information in your assay (incl. "non-variable" genes), is robust across experiments, learns signatures of transcriptional activity, and is colinear and non-negative (interpretable) rather than orthogonal and signed (not interpretable)

Singlet internally provides the **fastest implementation of NMF**. Cross-validation can take a few minutes for datasets with a few ten thousand cells, but is extremely scalable and runs excellently on HPC nodes and average laptops alike.
