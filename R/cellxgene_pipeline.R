#' Learn an NMF model from a cellxgene Seurat object
#'
#' @description Provide a link to download a cellxgene Seurat object, and this pipeline will return a standardized annotated NMF object at the optimal rank
#'
#' @details
#' This pipeline runs the following steps:
#' 1. Download a Seurat v4 object from the provided URL
#' 2. Preprocess the data and run NMF using parameters specified in the \code{...} argument
#' 3. Annotate the NMF model against existing multi-level factor information
#' 4. Extract the model and annotations and save to an RDS file
#'
#' @param url download url for a Seurat v4 object
#' @param ... arguments to \code{RunNMF}
#' @export
#' @md
#'
cellxgene_pipeline <- function(url, reps = 1, verbose = 3, L1 = 0.05, ...) {
    cat("reading ", url, "\n")
    A <- readRDS(url(url))
    if ("RNA" %in% names(A@assays)) {
        A@assays$RNA@key <- "RNA_"
        # keep only RNA assay
        A@assays <- list("RNA" = A@assays$RNA)
        cat(" normalizing...\n")
        A <- PreprocessData(A)
        cat(" running NMF...\n")
        t1 <- system.time({
            A <- RunNMF(A, reps = reps, verbose = 3, L1 = L1, ...)
        })[[3]]
        cat(" annotating NMF model...\n")
        A <- AnnotateNMF(A)

        model <- list(
            "w" = as(A@reductions$nmf@feature.loadings, "dgCMatrix"),
            "d" = A@reductions$nmf@stdev,
            "h" = as(A@reductions$nmf@cell.embeddings, "dgCMatrix"),
            "misc" = A@reductions$nmf@misc,
            "metadata" = A@meta.data,
            "dataset" = A@misc$title,
            "runtime" = t1
        )

        filename <- paste0(gsub("[^a-zA-Z]", "", A@misc$title), ".rds")
        cat(" saving model...\n")
        saveRDS(model, filename)
    }
}
