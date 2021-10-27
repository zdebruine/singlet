# singlet v.0.0.99

Singlet is an R toolkit for fast and simple single-cell analysis:
* clustering
* dimensional reduction

For vignettes, examples, and documentation, see https://zachdebruine.com/singlet

Singlet stable releases are publicly available on CRAN. Development versions on GitHub may be unstable. Singlet is being actively developed.

Get started: `install.packages("singlet")`

### What singlet can do:

cluster cells by alternating division/agglomeration from raw counts, find nearest neighbors by recursive bipartitioning from raw counts, and visualize on a UMAP reduction:
UMAP plot, heatmap of marker genes in clusters

Run non-negative matrix factorization on raw counts, find coordinated gene activities.

Transfer metadata from cells in one experiment onto cells in another experiment.

basic data QC/preprocessing

## Why singlet?
* fastest-in-class performance
* simple and intuitive algorithms
* simple functions returning simple results
* friendly to memory (your computers memory, that is)
* pretty plots

# Dependencies
* `RcppEigen` C++ library for all heavy lifting
* `ggplot2` for visualization
* `Matrix` for sparse matrix support

Absolutely everything else is done in-house in C++. YAY!
