# Small synthetic count data with `k` planted programs, so rank search has something to find.
make_counts <- function(genes = 300, cells = 120, k = 3, seed = 1) {
  set.seed(seed)
  W <- matrix(rexp(genes * k, 1 / 3) * (runif(genes * k) < 0.25), genes, k)
  H <- matrix(0, k, cells)
  grp <- rep(seq_len(k), length.out = cells)
  for (j in seq_len(cells)) H[grp[j], j] <- 1
  H <- H + matrix(runif(k * cells) * 0.1, k, cells)
  L <- W %*% H
  m <- matrix(rpois(length(L), L), genes, cells,
              dimnames = list(paste0("g", seq_len(genes)), paste0("c", seq_len(cells))))
  methods::as(m, "CsparseMatrix")
}

# A normalized Seurat object whose assay is a v3/v4 `Assay` or a v5 `Assay5`.
make_seurat <- function(assay_version = c("v5", "v3"), ...) {
  assay_version <- match.arg(assay_version)
  old <- options(Seurat.object.assay.version = assay_version)
  on.exit(options(old))
  obj <- SeuratObject::CreateSeuratObject(make_counts(...))
  Seurat::NormalizeData(obj, verbose = FALSE)
}
