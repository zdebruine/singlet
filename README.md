# singlet v.0.0.99

Singlet is an R toolkit for single-cell data analysis using non-negative matrix factorization.

## Dimensional Reduction with NMF

Analyze your single-cell assay with NMF:

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

![Integration with NMF](https://github.com/zdebruine/singlet/blob/main/readme_figures/Picture2.png)

Unlike Seurat anchor-based methods, integration with LNMF preserves unique signal and thus allows you to understand both the **shared and unique** signals in your different modalities/experiments.

Unlike LIGER integrative NMF, integration with LNMF perfectly separates unique signal and does not assume that shared and unique signals are of equal rank and linearly and additively correspond to one another.

LNMF falls short when a joint model does not capture significant overlap between groups. Work is ongoing to provide initializations to LNMF that address this.

## Ongoing Work

Singlet is being actively developed, thanks to funding from the Chan Zuckerberg Biohub:

* A new single-cell data class that uses 10x less memory than SCE or Seurat (and much faster)
* Full support for Seurat and SingleCellExperiment classes
* Out-of-core dimensional reduction with NMF
* Regularization and weighting to enable discovery of robust transcriptional signatures with NMF
* Spatially-aware dimensional reduction
* Extremely fast divisive clustering
