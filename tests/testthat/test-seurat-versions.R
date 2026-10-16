skip_if_not_installed("Seurat")
skip_if_not_installed("SeuratObject")

for (v in c("v3", "v5")) {
  local({
    version <- v
    test_that(paste0("RunNMF works on a Seurat ", version, " assay"), {
      obj <- make_seurat(version)
      expected <- if (version == "v5") "Assay5" else "Assay"
      expect_s4_class(obj[["RNA"]], expected)
      out <- RunNMF(obj, k = 3L, verbose = 0, seed = 1)
      expect_true("nmf" %in% names(out@reductions))
      emb <- SeuratObject::Embeddings(out, "nmf")
      expect_equal(dim(emb), c(ncol(obj), 3L))
      expect_true(all(emb >= 0))
    })

    test_that(paste0("RunNMF searches ranks, including rank 1, on a Seurat ", version, " assay"), {
      out <- RunNMF(make_seurat(version), k = 1:4, verbose = 0, seed = 1)
      expect_true("nmf" %in% names(out@reductions))
      expect_setequal(unique(out@reductions$nmf@misc$cv_data$k), 1:4)
    })

    test_that(paste0("RunNMF picks its own rank on a Seurat ", version, " assay"), {
      out <- RunNMF(make_seurat(version), verbose = 0, seed = 1)
      expect_gte(ncol(SeuratObject::Embeddings(out, "nmf")), 2L)
    })

    test_that(paste0("RunNMF can use the variable features of a Seurat ", version, " object"), {
      obj <- Seurat::FindVariableFeatures(make_seurat(version), nfeatures = 100, verbose = FALSE)
      out <- RunNMF(obj, k = 3L, features = "var.features", verbose = 0, seed = 1)
      expect_equal(nrow(SeuratObject::Loadings(out, "nmf")), 100L)
    })

    test_that(paste0("ProjectData projects onto a model on a Seurat ", version, " assay"), {
      obj <- make_seurat(version)
      fit <- RunNMF(obj, k = 3L, verbose = 0, seed = 1)
      w <- SeuratObject::Loadings(fit, "nmf")
      proj <- ProjectData(obj, w = w, reduction.name = "proj")
      expect_true("proj" %in% names(proj@reductions))
      expect_equal(dim(SeuratObject::Embeddings(proj, "proj")), c(ncol(obj), 3L))
    })

    test_that(paste0("an un-normalized Seurat ", version, " object is normalized on the way in"), {
      old <- options(Seurat.object.assay.version = version)
      on.exit(options(old), add = TRUE)
      raw <- SeuratObject::CreateSeuratObject(make_counts())
      out <- suppressWarnings(RunNMF(raw, k = 3L, verbose = 0, seed = 1))
      expect_true("nmf" %in% names(out@reductions))
    })

    test_that(paste0("PreprocessData log-normalizes a Seurat ", version, " object and honours scale.factor"), {
      old <- options(Seurat.object.assay.version = version)
      on.exit(options(old), add = TRUE)
      raw <- SeuratObject::CreateSeuratObject(make_counts())
      a <- singlet:::.assay_matrix(PreprocessData(raw, scale.factor = 1e4), "RNA", "data")
      b <- singlet:::.assay_matrix(PreprocessData(raw, scale.factor = 1e6), "RNA", "data")
      expect_s4_class(a, "dgCMatrix")
      expect_true(any(a@x != round(a@x)))        # no longer integer counts
      expect_gt(sum(b@x), sum(a@x))              # a larger scale factor gives larger values
      expect_equal(dim(a), dim(SeuratObject::GetAssayData(raw, assay = "RNA", layer = "counts")))
    })
  })
}

test_that("an assay with no usable layer says what to do instead of failing obscurely", {
  skip_if_not(utils::packageVersion("SeuratObject") >= "5.0.0")
  old <- options(Seurat.object.assay.version = "v5")
  on.exit(options(old), add = TRUE)
  obj <- SeuratObject::CreateSeuratObject(make_counts())
  obj[["RNA"]] <- split(obj[["RNA"]], f = rep(c("a", "b"), length.out = ncol(obj)))
  expect_error(suppressWarnings(RunNMF(obj, k = 3L, verbose = 0)), "JoinLayers|layer", ignore.case = TRUE)
})

test_that("the metadata method for Seurat objects no longer refers to an undefined variable", {
  skip_if_not_installed("SingleCellExperiment")
  obj <- make_seurat("v5")
  expect_no_error(S4Vectors::metadata(obj))
})
