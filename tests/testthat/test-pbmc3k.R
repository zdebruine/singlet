test_that("Testing pbmc3k data set",
          {
            data("pbmc3k", package="singlet")
            expect_true(all(c("i", "p", "Dim", "Dimnames", "x", "cell_type")
                            %in% names(pbmc3k)))
            expect_equal(length(pbmc3k$cell_type), pbmc3k$Dim[[2]])
          })
