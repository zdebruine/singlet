# singlet v.0.0.99

Singlet is an R toolkit for single-cell data analysis using non-negative matrix factorization.

## Dimensional Reduction with NMF

NMF can be used to do almost anything that PCA can do. In addition, NMF:
* imputes missing signal
* always has an optimal rank
* uses all the information in your assay (incl. "non-variable" genes)
* is robust across experiments
* learns signatures of transcriptional activity
* is colinear, rather than orthogonal

Singlet directly provides the absolute fastest implementation of NMF. Cross-validation can take a few minutes for datasets with a few ten thousand cells, but is extremely scalable and runs excellently on HPC nodes or average laptops alike.

## Integration with Linked NMF

## Ongoing Work

Singlet is being actively developed, thanks to funding from the Chan Zuckerberg Biohub:

* A new single-cell data class that uses 10x less memory than SCE or Seurat (and much faster)
* Full support for Seurat and SingleCellExperiment classes
* Out-of-core dimensional reduction with NMF
* Regularization and weighting to enable discovery of robust transcriptional signatures with NMF
* Spatially-aware dimensional reduction
* Extremely fast divisive clustering
