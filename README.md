# singlet v.0.0.99

Singlet is an R toolkit for fast, robust, and interpretable single-cell experiment analysis:

* dimensional reduction with NMF
* dataset integration
* multimodal integration
* biological process discovery
* clustering (divisive and graph-based clustering)
* spatial RNA deconvolution

## Plan

End-to-end single-cell analysis package capable of scaling to the largest datasets being generated today.

Performant sparse matrix operations, dimensional reduction, and clustering. R S4 front-end with backend largely in Rcpp, and Eigen C++ BLAS for linear algebra operations.

Likely not a Bioconductor package, because the Bioconductor S4 ecosystem is very unfriendly to C++ wrappers.


Run non-negative matrix factorization on raw counts, find coordinated gene activities.

## Development timeline

Before Bioc2022:  Package skeleton with automated NMF dimensional reduction, UMAP visualization, graph-based and divisive clustering, factor annotation against GO terms and any sample metadata in the assay.

Late 2022: Add support for multi-experiment integration (identify and remove batch effect factors) and multi-modal integration.

Publish first major manuscript for the project.

2023: Meta-analysis of existing data across many experiments and modalities to learn master models of biological information. Use this information to annotate single-cell RNA, CITE, ATAC assays.
