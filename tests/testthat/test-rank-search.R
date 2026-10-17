test_that("cross_validate_nmf accepts rank 1 (a one-row index must not drop to a vector)", {
  A <- make_counts()
  cv <- cross_validate_nmf(A, ranks = 1:3, n_replicates = 1, verbose = 0)
  expect_s3_class(cv, "cross_validate_nmf_data")
  expect_setequal(unique(cv$k), 1:3)
})

test_that("cross_validate_nmf works on a dense matrix too", {
  cv <- cross_validate_nmf(as.matrix(make_counts(100, 60)), ranks = 1:2, n_replicates = 1, verbose = 0)
  expect_setequal(unique(cv$k), 1:2)
})

test_that("ard_nmf finds a rank and returns a model of that rank", {
  m <- ard_nmf(make_counts(), verbose = 0)
  expect_gte(ncol(m$w), 2L)
  expect_equal(ncol(m$w), nrow(m$h))
  expect_true(all(m$w >= 0) && all(m$h >= 0))
  expect_s3_class(m$cv_data, "cross_validate_nmf_data")
})

test_that("a matrix with no signal does not crash the loop with 'missing value where TRUE/FALSE needed'", {
  z <- methods::as(matrix(0, 40, 30, dimnames = list(paste0("g", 1:40), paste0("c", 1:30))), "CsparseMatrix")
  res <- suppressWarnings(tryCatch(ard_nmf(z, verbose = 0), error = function(e) e))
  if (inherits(res, "error")) {
    expect_match(conditionMessage(res), "rank search failed|finite", ignore.case = TRUE)
    expect_false(grepl("TRUE/FALSE|length zero", conditionMessage(res)))
  }
  expect_warning(tryCatch(ard_nmf(z, verbose = 0), error = function(e) NULL), "not finite")
})

test_that("GetBestRank on nothing fitted is NA, not an error", {
  empty <- data.frame(k = integer(), rep = integer(), test_error = double(), iter = integer())
  expect_true(is.na(GetBestRank(empty)))
})

test_that("a tiny sparse matrix still gets a rank", {
  set.seed(2)
  tiny <- methods::as(matrix(rpois(60 * 20, 0.3), 60, 20,
                             dimnames = list(paste0("g", 1:60), paste0("c", 1:20))), "CsparseMatrix")
  expect_gte(ncol(ard_nmf(tiny, verbose = 0)$w), 2L)
})
