# singlet v.0.0.99

**singlet is being actively developed (including documentation) and there are still bugs in the LNMF function**.

Singlet brings fast Non-negative Matrix Factorization, new integration methods, and more to every Seurat user.

## Get Started!

```{R}
install.packages("singlet")
library(singlet)
```

Planned vignettes:
* Guided clustering tutorial with `SeuratData::pbmc3k` dataset
* Batch integration with `SeuratData::ifnb` dataset (PBMC +/- stimulation)
* Multi-modal integration with `SeuratData::bmcite`

## Dimensional Reduction with NMF

Analyze your single-cell assay with NMF:

```{R}
library(singlet)
library(Seurat)
library(SeuratData)
library(dplyr)
library(cowplot)
data(pbmc3k)
pbmc3k <- NormalizeData(pbmc3k)
set.seed(123)
pbmc3k <- RunNMF(pbmc3k, k = seq(2, 20, 2), n_replicates = 1)
pbmc3k <- RunUMAP(pbmc3k, reduction = "nmf", dims = 1:ncol(pbmc3k@reductions$nmf))
plot_grid(RankPlot(pbmc3k) + NoLegend(), DimPlot(pbmc3k) + NoLegend(), labels = "auto", ncol = 2)
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

## Integration with Linked NMF

Learn an integrated model of information across modalities or sample batches and explore **shared and unique** signals in each of your groups.

```{R}
library(singlet)
library(SeuratData)
data(ifnb)
ifnb <- NormalizeData(ifnb)
ifnb <- RunNMF(ifnb, k = 30, split.by = "stim")
ifnb <- RunLNMF(ifnb, split.by = "stim")
jnmf_plot <- MetadataPlot(ifnb, split.by = "stim", reduction = "nmf")
lnmf_plot <- MetadataPlot(ifnb, split.by = "stim", reduction = "lnmf")
ifnb <- RunUMAP(ifnb, reduction = "nmf", dims = 1:ncol(ifnb@reductions$nmf), reduction.name = "jnmf_all")
ifnb <- RunUMAP(ifnb, reduction = "lnmf", dims = 1:ncol(ifnb@reductions$lnmf), reduction.name = "lnmf_all")
ifnb <- RunUMAP(ifnb, reduction = "lnmf", dims = GetSharedFactors(ifnb, split.by = "stim"), reduction.name = "lnmf_shared")
selected_factors <-  which(rowSums(ifnb@reductions$lnmf@misc$link_matrix == 0) == 0)
ifnb <- RunUMAP(ifnb, reduction = "nmf", dims = selected_factors, reduction.name = "jnmf_shared")
plot_grid(jnmf_plot + ggtitle("joint NMF factor/batch weights") + theme(legend.position = "none"), lnmf_plot + ggtitle("linked NMF factor/batch weights") + theme(legend.position = "none"), get_legend(jnmf_plot), rel_widths = c(1, 1, 0.4), ncol = 3)

p_jnmf_umap <- DimPlot(ifnb, reduction = "jnmf_all", group.by = "stim")
p_jnmf_shared_umap <- DimPlot(ifnb, reduction = "jnmf_shared", group.by = "stim")
p_lnmf_all_umap <- DimPlot(ifnb, reduction = "lnmf_all", group.by = "stim")
p_lnmf_shared_umap <- DimPlot(ifnb, reduction = "lnmf_shared", group.by = "stim")

plot_grid(plot_grid(
  p_jnmf_umap + ggtitle("joint NMF") + theme(legend.position = "none") + theme(legend.position = "none"), 
  p_jnmf_shared_umap + ggtitle("joint NMF, selected factors") + theme(legend.position = "none"), 
  p_lnmf_all_umap + ggtitle("linked NMF, all factors") + theme(legend.position = "none"), 
  p_lnmf_shared_umap + ggtitle("linked NMF, shared factors") + theme(legend.position = "none"), 
  nrow = 2, ncol = 2
 ), get_legend(p_jnmf_umap), ncol = 2, rel_widths = c(1, 0.2))
```

Because unique and shared features are now completely separated, we can run GSEA on these features to better understand whether they are technical "batch effects" or biological differences:

```
ifnb <- RunGSEA(ifnb, reduction = "lnmf", dims = GetUniqueFactors(ifnb, split.by = "stim"))
GSEAHeatmap(ifnb, reduction = "lnmf")
```

![Integration with NMF](https://github.com/zdebruine/singlet/blob/main/readme_figures/Picture2.png)

Unlike Seurat anchor-based methods, integration with LNMF preserves unique signals in the reduction and thus allows you to understand both the **shared and unique** signals in your different modalities/experiments.

Unlike LIGER integrative NMF, integration with LNMF completely uncouples shared and unique signals from one another and does not assume that both shared and unique signals in any given factor can be mapped to a given sample by the same weight (and thus correspond linearly and additively to one another).

Why not LNMF? If signals are highly disparate, they will not be aligned. Work is ongoing to provide initializations to LNMF that can align even the most disparate signals. However, care should be taken to properly interpret results and avoid "forcing" too much integration when signals are radically different.

## Ongoing Work

Singlet is being actively developed:

* A new single-cell data class that uses 10x less memory than SCE or Seurat (and much faster)
* Full support for Seurat and SingleCellExperiment classes
* Out-of-core dimensional reduction with NMF
* Regularization and weighting to enable discovery of robust transcriptional signatures with NMF
* Spatially-aware dimensional reduction
* Graph-regularization for more aligned joint NMF models to initialize LNMF
