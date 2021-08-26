# singlet

Easy and fast single-cell analysis in R for decomposing coordinated gene activities and classifying cell identities, without devils in the details.

Get started: `install.packages("singlet")`

### What singlet can do:

cluster cells by alternating division/agglomeration from raw counts, find nearest neighbors by recursive bipartitioning from raw counts, and visualize on a UMAP reduction:
UMAP plot, heatmap of marker genes in clusters

Run non-negative matrix factorization on raw counts, find coordinated gene activities, map paracrine signaling activity between factors, and plot paracrine signaling networks on a UMAP reduction:
PLOTS

Transfer metadata from cells in one experiment onto cells in another experiment.

basic data QC/preprocessing

## Why singlet?
* very fast. Almost no patience needed
* algos are simple, almost no hyperparameters
* simple functions return simple results
* friendly to memory
* pretty plots
* normalization, scaling, and centering are NEVER needed. Raw counts to the moon!
* always uses all genes, because only highly variable genes are NOT ENOUGH

# Dependencies
* `RcppEigen` C++ library for all heavy lifting
* `ggplot2` for visualization
* `Matrix` for sparse matrix support

Absolutely everything else is done in-house in C++. YAY!
