test_that("Testing pbmc3k data set",
          {
            data("pbmc3k", package="singlet")
            all(c("i", "p", "Dim", "Dimnames", "x", "cell_type") 
                %in% names(pbmc3k))
            expect_true(TRUE)
          })
