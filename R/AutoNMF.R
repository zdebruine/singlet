#' Automatic Non-negative Matrix Factorization
#'
#' AutoNMF quickly learns an NMF model on a matrix at the largest rank for which overfitting does not occur, which usually is the rank that maximizes imputation accuracy and generalization of the model to other datasets.
#' In other words, AutoNMF learns a maximally robust and minimally biased model.
#'
#' @details
#' AutoNMF uses a combination of downsampling, coordinate descent, and detection of overfitting to quickly determine the optimal rank of an NMF model. This automatic rank determination procedure is the following:
#'  1. Select a random sample of columns, where the sample size is a power of 2 given by \code{auto.subsample.exponent = 10} by default.  For example, \code{2^10} corresponds to 1024 samples.
#'  2. Select an initial rank (default \code{auto.k.init = 2}) and a step size (default \code{auto.step.size = 2}). Using this algorithm, it will be impossible to decrease the rank from the initial value provided.
#'  3. Initialize \code{w} with values in a random uniform distribution.  Define a test set as a random speckled pattern with some density (default \code{auto.test.set.density = 0.05}).
#'  4. Fit the model at the initial rank, masking test set values. Stop alternating least squares updates when any of three stopping criterion are satisfied:
#'      1) a maximum number of iterations (default \code{auto.maxit = 100})
#'      2) Pearson correlation distance between models across consecutive iterations (default \code{auto.tol = 1e-4})
#'      3) significant increase in test set reconstruction error over the minimum test set reconstruction error during model fitting (default \code{auto.overfit.tol = 1e-3})
#'  5. If the model fitting was not terminated due to overfitting, record the test set reconstruction error at convergence. Increase the rank by \code{auto.step.size}, refit, and repeat this procedure until overfitting is detected.
#'  6. If the model is overfit, increase the column sample size by a power of 2 and refit.
#'  7. Once the sample size cannot be increased by a power of 2 without exceeding the number of samples in the dataset, reduce the rank by \code{auto.step.size} from the last detected overfit rank, remove the test set, and factorize the full dataset.
#'
#' The seed may be specified directly to the function using \code{use.seed}. If \code{use.seed = NULL}, the seed will be derived from the current R environment \code{.Random.seed}. The random seed will control downsampling, test set membership, and initialization of the model.
#'
#' @md
#' @details
#' @import dplyr
#' @param object a dense or sparse matrix to factorize
#' @param reps number of times to repeat the automatic rank determination procedure before learning the final model once on the consensus best rank
#' @param tol stopping criteria for the final model fit. Given as the Pearson correlation distance between \code{w} across consecutive iterations.
#' @param maxit stopping criteria for the final model fit. Maximum number of alternating least squares iterations.
#' @param L1 L1/LASSO regularization penalty in the range (0, 1] to increase the sparsity of each factor in the model.
#' @param L2 L2/Ridge regularization penalty in the range (0, 1] to increase the angle between factors in the model.
#' @param threads number of threads to use.
#' @param ... advanced parameters with reasonable defaults that should not normally need to be changed:
#' \code{auto.subsample.exponent = 10},
#' \code{auto.k.init = 2},
#' \code{auto.step.size = 2},
#' \code{auto.test.set.density = 0.05},
#' \code{auto.maxit = 100},
#' \code{auto.overfit.tol = 1e-3},
#' \code{auto.tol = 1e-4}
#'
#' @export
#'
AutoNMF <- function(object, reps = 3, tol = 1e-5, maxit = 100, verbose = 2, L1 = 0.01, L2 = 0, threads = 0, use.seed = NULL, ...) {
    stopifnot("L1 penalty must be strictly in the range (0, 1]" = L1 < 1 & L1 >= 0)
    stopifnot("L2 penalty must be strictly in the range (0, 1]" = L2 < 1 & L2 >= 0)

    # set advanced fit parameters
    dot_params <- list(...)
    if (is.null(dot_params$auto.step.size)) dot_params$auto.step.size <- 2
    if (is.null(dot_params$auto.test.set.density)) dot_params$auto.test.set.density <- 0.05
    inv_test_set_density <- ceiling(1 / dot_params$auto.test.set.density)
    if (is.null(dot_params$auto.maxit)) dot_params$auto.maxit <- 100
    if (is.null(dot_params$auto.overfit.tol)) dot_params$auto.overfit.tol <- 1e-3
    if (is.null(dot_params$auto.tol)) dot_params$auto.tol <- 1e-4

    # set random seed
    if (is.null(use.seed)) {
        use.seed <- abs(.Random.seed[[3]])
        set.seed(random_seed)
    }

    object_t <- t(object)

    # calculate downsample vector
    max_exponent <- floor(log2(ncol(object)))
    auto_subsample_exponent <- ifelse(is.null(dot_params$auto.subsample.exponent), 10, dot_params$auto.subsample.exponent)
    if(max_exponent < auto_subsample_exponent) max_exponent <- auto_subsample_exponent
    subsample_size <- 2^(auto_subsample_exponent:max_exponent)
    subsample_size[[length(subsample_size)]] <- ncol(object)

    # Automatic rank determination procedure
    cv_data <- list()
    best_ranks <- c()
    for (current_rep in 1:reps) {
        use.seed <- use.seed + 1
        if (verbose >= 1 && n_replicates > 1) cat("\nREPLICATE", curr_rep, "/", n_replicates, "\n")
        current_rank <- ifelse(is.null(dot_params$auto.k.init), 0, dot_params$auto.k.init - dot_params$auto.step.size)
        for (sample_size in subsample_size) {
            is_overfit <- FALSE
            while (!is_overfit) {
                current_rank <- current_rank + dot_params$auto.step.size
                if (verbose > 0) cat("sample size:", sample_size, ", k:", current_rank, ", replicate:", current_rep, "\n")
                # TO DO: write subsample ard_nmf using iterators
                model <- c_ard_nmf_subsample(object, object_t, n_samples, current_rank, dot_params$auto.tol, dot_params$auto.maxit, verbose > 2, L1, L2, threads, use.seed, inv_test_set_density, dot_params$auto.overfit.tol, 1)
                is_overfit <- tail(model$overfit_score, 1L) > dot_params$auto.overfit.tol
                cv_data[[length(cv_data) + 1]] <- data.frame(
                    "rep" = current_rep,
                    "sample_size" = sample_size,
                    "k" = current_rank,
                    "test_error" = model$test_mse,
                    "it" = model$iter,
                    "tol" = model$tol,
                    "overfit_tol" = overfit_score,
                    "is_overfit" = is_overfit
                )
            }
        }
        best_ranks <- c(best_ranks, current_rank)
    }
    cv_data <- do.call(rbind, cv_data)
    best_rank <- floor(mean(best_ranks))

    # learn final nmf model
    if (verbose > 0) cat("\nFitting final model at k =", best_rank, "\n")
    w_init <- r_matrix(current_rank, nrow(object), use.seed)
    model <- c_nmf(object, object_t, tol, maxit, verbose > 2, L1, L2, threads, w_init)

    # format and return model
    if (transposed) {
        w <- model$w
        model$w <- model$h
        model$h <- w
    }
    model$cv_data <- cv_data
    sort_index <- order(model$d, decreasing = TRUE)
    model$d <- model$d[sort_index]
    model$w <- t(model$w)[, sort_index]
    model$h <- model$h[sort_index, ]
    rownames(model$w) <- rownames(object)
    colnames(model$h) <- colnames(object)
    colnames(model$w) <- rownames(model$h) <- paste0("NMF_", 1:ncol(model$w))
    model
}
