# Seurat object access that works across versions.
#
# Seurat v3/v4 keep a normalized matrix in `Assay@data`; v5 keeps named layers in
# `Assay5`, and `@data` / `$data` / `@var.features` are not there or not stable.
# Everything in this package that reads or writes an assay matrix goes through
# the helpers below, which use the SeuratObject layer API (`layer =` on
# SeuratObject >= 5, `slot =` before it).

.so5 <- function() utils::packageVersion("SeuratObject") >= "5.0.0"

# The matrix in `layer` ("data" or "counts") of `assay`, as a dgCMatrix, or NULL when
# the assay has no such layer (or its layers are split and must be joined first).
.assay_matrix <- function(object, assay, layer = "data") {
  m <- tryCatch(
    if (.so5()) {
      SeuratObject::GetAssayData(object, assay = assay, layer = layer)
    } else {
      SeuratObject::GetAssayData(object, assay = assay, slot = layer)
    },
    error = function(e) NULL
  )
  if (is.null(m) || is.null(dim(m)) || any(dim(m) == 0L)) return(NULL)
  if (!inherits(m, "dgCMatrix")) m <- as(as(as(m, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  m
}

.set_assay_matrix <- function(object, assay, layer, value) {
  if (.so5()) {
    SeuratObject::SetAssayData(object, assay = assay, layer = layer, new.data = value)
  } else {
    SeuratObject::SetAssayData(object, assay = assay, slot = layer, new.data = value)
  }
}

# Returns `object` with log-normalized data in the assay's data layer: leaves it alone when
# it already is normalized, otherwise normalizes the counts (an integer-valued "data" layer
# means the object was never normalized).
.ensure_normalized <- function(object, assay) {
  d <- .assay_matrix(object, assay, "data")
  if (!is.null(d) && length(d@x) > 0L && sum(as.integer(d@x)) != sum(d@x)) return(object)
  PreprocessData(object, assay = assay)
}
