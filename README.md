# singlet v.0.0.99

Singlet brings fast Non-negative Matrix Factorization, new integration methods, and more to every Seurat user.

Now you can run NMF with Seurat. Why still use PCA?

## Get Started!

```{R}
install.packages("singlet")
library(singlet)
```

Vignettes:
* Guided clustering tutorial with `SeuratData::pbmc3k` dataset
* Batch integration with `SeuratData::ifnb` dataset (PBMC +/- stimulation)
* Multi-modal integration with `SeuratData::bmcite`

## Dimensional Reduction with NMF

Analyze your single-cell assay with NMF:

```{R}
library(singlet)
library(SeuratData)
data(pbmc3k)
pbmc3k <- pbmc3k %>% 
  LogNormalize() %>% 
  RunNMF(k = 2:30, n_replicates = 3) %>% 
  RunUMAP()
cowplot::plot_grid(RankPlot(pbmc3k), DimPlot(pbmc3k), labels = "auto", ncol = 3)
```

![NMF workflow](https://github.com/zdebruine/singlet/blob/main/readme_figures/Picture1.png)

NMF can do almost anything that PCA can do, but also:
* imputes missing signal
* always has an optimal rank (for variance-stabilized data)
* uses all the information in your assay (incl. "non-variable" genes)
* is robust across experiments
* learns signatures of transcriptional activity
* is colinear (interpretable), rather than orthogonal (not interpretable)

Singlet directly provides the **absolute fastest implementation of NMF**. Cross-validation can take a few minutes for datasets with a few ten thousand cells, but is extremely scalable and runs excellently on HPC nodes and average laptops alike.

## Integration with Linked NMF

Learn an integrated model of information across modalities or experiments and explore **shared and unique** signals in each of your groups.

```{R}
library(singlet)
library(SeuratData)
data(ifnb)
ifnb <- ifnb %>% 
  singlet::NormalizeData() %>% 
  singlet::RunNMF(k = 30, split.by = "stim") %>% 
  singlet::RunLNMF(split.by = "stim")
MetadataPlot(ifnb, split.by = "stim")
ifnb <- RunUMAP(ifnb, reduction = "lnmf", dims = GetSharedFactors(ifnb))
DimPlot(ifnb)
```

![Integration with NMF](https://github.com/zdebruine/singlet/blob/main/readme_figures/Picture2.png)

Unlike Seurat anchor-based methods, integration with LNMF preserves unique signals in the reduction and thus allows you to understand both the **shared and unique** signals in your different modalities/experiments.

Unlike LIGER integrative NMF, integration with LNMF completely "unties" shared and unique signals from one another and does not assume that both shared and unique signals in any given factor can be mapped to a given sample by the same weight (and thus correspond linearly and additively to one another).

Why not LNMF? If signals are highly disparate, they will not be aligned. Work is ongoing to provide initializations to LNMF that can align even the most disparate signals. However, care should be taken to properly interpret results and avoid "forcing" too much integration when signals are radically different.

## Ongoing Work

Singlet is being actively developed:

* A new single-cell data class that uses 10x less memory than SCE or Seurat (and much faster)
* Full support for Seurat and SingleCellExperiment classes
* Out-of-core dimensional reduction with NMF
* Regularization and weighting to enable discovery of robust transcriptional signatures with NMF
* Spatially-aware dimensional reduction
* Extremely fast divisive clustering
