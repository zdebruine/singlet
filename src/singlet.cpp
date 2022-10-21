#define EIGEN_NO_DEBUG
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

//[[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

#include <singlet.h>

// RANDOM NUMBER GENERATOR  ---------------------------------------------------------------------------------
// two-dimensional linear congruential random number generator using
//   Marsaglia's xorshift32 algorithm together with xoroshiro128++ and a
//   random seed. The generated RNG is non-correlated with sequences in `i` and `j`.
//   The returned logical is a probability of 1/max and is transpose-identical (independent of order of `i` and `j`).
//   "t" implicitly transposes the matrix
inline bool is_masked(uint32_t i, uint32_t j, uint32_t state, uint32_t max, const bool t) {
    if (t) std::swap(i, j);
    // if (j >= i) std::swap(i, j);
    uint64_t ij = (i + 1) * (i + 2) / 2 + j + 1;
    ij ^= ij << 13 | (i << 17);
    ij ^= ij >> 7 | (j << 5);
    ij ^= ij << 17;
    uint64_t s = state ^ ij;
    s ^= s << 23;
    s = s ^ ij ^ (s >> 18) ^ (ij >> 5);
    s += ij;
    return s % max == 0;
}

// MATRIX NORMALIZATION BY GROUP -----------------------------------------------------------------------------
// simply makes all values in each group of samples sum to the same value
//[[Rcpp::export]]
Rcpp::S4 weight_by_split(const Rcpp::S4& A_, Rcpp::IntegerVector split_by, const size_t n_groups) {
    Rcpp::SparseMatrix A(A_);
    A.clone();

    // get sum of each group
    Rcpp::NumericVector sums(n_groups);
    for (size_t j = 0; j < split_by.size(); ++j) {
        for (Rcpp::SparseMatrix::InnerIterator it(A, j); it; ++it) {
            sums[split_by[j]] += it.value();
        }
    }

    // calculate multiplication factor for each group
    for (size_t j = 1; j < sums.size(); ++j)
        sums[j] /= sums[0];

    // normalize each column using the group sum relative to the first group sum
    for (size_t i = 0; i < split_by.size(); ++i) {
        if (split_by[i] != 0) {
            for (Rcpp::SparseMatrix::InnerIterator it(A, i); it; ++it)
                it.value() /= sums[split_by[i]];
        }
    }

    return A.wrap();
}

//[[Rcpp::export]]
Eigen::MatrixXd rowwise_compress(Rcpp::SparseMatrix& A, const size_t n = 10, const size_t threads = 0) {
    const size_t n_rows = (size_t)std::floor(A.rows() / n);
    Eigen::MatrixXd res(n_rows, A.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t col = 0; col < A.cols(); ++col)
        for (Rcpp::SparseMatrix::InnerIterator it(A, col); it; ++it) {
            size_t row = (size_t)std::floor(it.row() / n);
            res(row, col) += it.value();
        }
    res.array() /= n;  // calculate mean
    return res;
}

// NMF HELPER FUNCTIONS ---------------------------------------------------------------------------------
// Pearson correlation between two matrices
inline double cor(Eigen::MatrixXd& x, Eigen::MatrixXd& y) {
    double x_i, y_i, sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
    const size_t n = x.size();
    for (size_t i = 0; i < n; ++i) {
        x_i = (*(x.data() + i));
        y_i = (*(y.data() + i));
        sum_x += x_i;
        sum_y += y_i;
        sum_xy += x_i * y_i;
        sum_x2 += x_i * x_i;
        sum_y2 += y_i * y_i;
    }
    return 1 - (n * sum_xy - sum_x * sum_y) / std::sqrt((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y));
}

// fast symmetric matrix multiplication, A * A.transpose()
inline Eigen::MatrixXd AAt(const Eigen::MatrixXd& A) {
    Eigen::MatrixXd AAt = Eigen::MatrixXd::Zero(A.rows(), A.rows());
    AAt.selfadjointView<Eigen::Lower>().rankUpdate(A);
    AAt.triangularView<Eigen::Upper>() = AAt.transpose();
    AAt.diagonal().array() += 1e-15;
    return AAt;
}

// subset columns from a matrix (deep copy)
// could not figure out how to do better than a deep copy:
//   https://stackoverflow.com/questions/72100483/matrix-multiplication-of-an-eigen-matrix-for-a-subset-of-columns
inline Eigen::MatrixXd submat(const Eigen::MatrixXd& x, const std::vector<uint32_t>& col_indices) {
    Eigen::MatrixXd x_(x.rows(), col_indices.size());
    for (size_t i = 0; i < col_indices.size(); ++i)
        x_.col(i) = x.col(col_indices[i]);
    return x_;
}

// scale rows in w (or h) to sum to 1 and put previous rowsums in d
void scale(Eigen::MatrixXd& w, Eigen::VectorXd& d) {
    d = w.rowwise().sum();
    d.array() += 1e-15;
    for (size_t i = 0; i < w.rows(); ++i)
        for (size_t j = 0; j < w.cols(); ++j)
            w(i, j) /= d(i);
};

// NNLS SOLVER OF THE FORM ax=b ---------------------------------------------------------------------------------
// optimized and modified from github.com/linxihui/NNLM "c_nnls" function
inline void nnls(Eigen::MatrixXd& a, Eigen::VectorXd& b, Eigen::MatrixXd& x, const size_t col) {
    double tol = 1;
    for (uint8_t it = 0; it < 100 && (tol / b.size()) > 1e-8; ++it) {
        tol = 0;
        for (size_t i = 0; i < x.rows(); ++i) {
            double diff = b(i) / a(i, i);
            if (-diff > x(i, col)) {
                if (x(i, col) != 0) {
                    b -= a.col(i) * -x(i, col);
                    tol = 1;
                    x(i, col) = 0;
                }
            } else if (diff != 0) {
                x(i, col) += diff;
                b -= a.col(i) * diff;
                tol += std::abs(diff / (x(i, col) + 1e-15));
            }
        }
    }
}

// NMF PROJECTION FUNCTIONS  ------------------------------------------------------------------------------

// update h given A and w
inline void predict(Rcpp::SparseMatrix A, const Eigen::MatrixXd& w, Eigen::MatrixXd& h, const double L1, const double L2, const int threads) {
    Eigen::MatrixXd a = AAt(w);
    if (L2 != 0) a.diagonal().array() *= (1 - L2);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < h.cols(); ++i) {
        if (A.p[i] == A.p[i + 1]) continue;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(h.rows());
        for (Rcpp::SparseMatrix::InnerIterator it(A, i); it; ++it)
            b += it.value() * w.col(it.row());
        b.array() -= L1;
        nnls(a, b, h, i);
    }
}

// update h given A and w
inline void predict(Eigen::MatrixXd A, const Eigen::MatrixXd& w, Eigen::MatrixXd& h, const double L1, const double L2, const int threads) {
    Eigen::MatrixXd a = AAt(w);
    if (L2 != 0) a.diagonal().array() *= (1 - L2);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < h.cols(); ++i) {
        Eigen::VectorXd b = A.col(i) * w;
        b.array() -= L1;
        nnls(a, b, h, i);
    }
}

//[[Rcpp::export]]
Rcpp::List c_project_model(Rcpp::SparseMatrix A, Eigen::MatrixXd w, const double L1, const double L2, const int threads) {
    if (w.rows() == A.rows()) w = w.transpose();
    Eigen::VectorXd d(w.rows());
    scale(w, d);
    Eigen::MatrixXd h(w.rows(), A.cols());
    predict(A, w, h, L1, L2, threads);
    scale(h, d);
    return (Rcpp::List::create(Rcpp::Named("h") = h, Rcpp::Named("d") = d));
}

// update h given A and w while masking values in h
inline void predict_link(Rcpp::SparseMatrix A, const Eigen::MatrixXd& w, Eigen::MatrixXd& h, const double L1, const double L2, const int threads,
                         const Eigen::MatrixXd& link_h) {
    Eigen::MatrixXd a = AAt(w);
    if (L2 != 0) a.diagonal().array() *= (1 - L2);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < h.cols(); ++i) {
        if (A.p[i] == A.p[i + 1]) continue;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(h.rows());
        for (Rcpp::SparseMatrix::InnerIterator it(A, i); it; ++it)
            b += it.value() * w.col(it.row());
        b.array() -= L1;
        for (size_t j = 0; j < link_h.rows(); ++j)
            b[j] *= link_h(j, i);
        nnls(a, b, h, i);
    }
}

// update h given A and w while masking a random speckled test set
inline void predict_mask(Rcpp::SparseMatrix A, const uint32_t seed, const uint32_t inv_density, const Eigen::MatrixXd& w,
                         Eigen::MatrixXd& h, const double L1, const double L2, const int threads, const bool mask_t) {
    Eigen::MatrixXd a = AAt(w);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < h.cols(); ++i) {
        if (A.p[i] == A.p[i + 1]) continue;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(h.rows());
        Rcpp::SparseMatrix::InnerIterator it(A, i);
        std::vector<uint32_t> idx;
        idx.reserve(A.rows() / inv_density);
        for (uint32_t j = 0; j < A.rows(); ++j) {
            if (is_masked(i, j, seed, inv_density, mask_t)) {
                idx.push_back(j);
                if (it && j == it.row())
                    ++it;
            } else if (it && j == it.row()) {
                b += it.value() * w.col(j);
                ++it;
            }
        }
        b.array() -= L1;
        Eigen::MatrixXd wsub = submat(w, idx);
        Eigen::MatrixXd asub = AAt(wsub);
        Eigen::MatrixXd a_i = a - asub;
        if (L2 != 0) a_i.diagonal().array() *= (1 - L2);
        nnls(a_i, b, h, i);
    }
}

// update h given A and w while masking a random speckled test set
inline void predict_mask(const Eigen::MatrixXd& A, const uint32_t seed, const uint32_t inv_density, const Eigen::MatrixXd& w,
                         Eigen::MatrixXd& h, const double L1, const double L2, const int threads, const bool mask_t) {
    Eigen::MatrixXd a = AAt(w);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < h.cols(); ++i) {
        Eigen::VectorXd b = Eigen::VectorXd::Zero(h.rows());
        std::vector<uint32_t> idx;
        idx.reserve(A.rows() / inv_density);
        for (uint32_t j = 0; j < A.rows(); ++j) {
            if (is_masked(i, j, seed, inv_density, mask_t)) {
                idx.push_back(j);
            } else {
                b += A(j, i) * w.col(j);
            }
        }
        b.array() -= L1;
        Eigen::MatrixXd wsub = submat(w, idx);
        Eigen::MatrixXd asub = AAt(wsub);
        Eigen::MatrixXd a_i = a - asub;
        if (L2 != 0) a_i.diagonal().array() *= (1 - L2);
        nnls(a_i, b, h, i);
    }
}

// CALCULATE ERROR OF TEST SET ----------------------------------------------------------------------------

// calculate mean squared error of the model at test set indices only
inline double mse_test(Rcpp::SparseMatrix A, const Eigen::MatrixXd& w, Eigen::VectorXd& d, Eigen::MatrixXd& h,
                       const uint32_t seed, const uint32_t inv_density, const uint16_t threads) {
    // multiply w by d
    Eigen::MatrixXd w_ = w.transpose();
    for (size_t i = 0; i < w.cols(); ++i)
        for (size_t j = 0; j < w.rows(); ++j)
            w_(i, j) *= d(j);

    Eigen::VectorXd losses(h.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (uint32_t j = 0; j < h.cols(); ++j) {
        uint32_t n = 0;
        double s = 0;
        Rcpp::SparseMatrix::InnerIterator it(A, j);
        for (uint32_t i = 0; i < A.rows(); ++i) {
            if (is_masked(j, i, seed, inv_density, false)) {
                ++n;
                if (it && i == it.row()) {
                    s += std::pow(w_.row(i) * h.col(j) - it.value(), 2);
                    ++it;
                } else {
                    s += std::pow(w_.row(i) * h.col(j), 2);
                }
            } else if (it && i == it.row()) {
                ++it;
            }
        }
        losses(j) = (n > 0) ? s / n : 0;
    }
    return losses.sum() / h.cols();
}

// calculate mean squared error of the model at test set indices only
inline double mse_test(const Eigen::MatrixXd& A, const Eigen::MatrixXd& w, Eigen::VectorXd& d, Eigen::MatrixXd& h,
                       const uint32_t seed, const uint32_t inv_density, const uint16_t threads) {
    // multiply w by d
    Eigen::MatrixXd w_ = w.transpose();
    for (size_t i = 0; i < w.cols(); ++i)
        for (size_t j = 0; j < w.rows(); ++j)
            w_(i, j) *= d(j);

    Eigen::VectorXd losses(h.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (uint32_t j = 0; j < h.cols(); ++j) {
        uint32_t n = 0;
        double s = 0;
        for (uint32_t i = 0; i < A.rows(); ++i) {
            if (is_masked(j, i, seed, inv_density, false)) {
                ++n;
                s += std::pow(w_.row(i) * h.col(j) - A(i, j), 2);
            }
        }
        losses(j) = (n > 0) ? s / n : 0;
    }
    return losses.sum() / h.cols();
}

// NMF FUNCTIONS ---------------------------------------------------------------------------------------
// no linking or masking
template <class Matrix>
Rcpp::List c_nmf_base(Matrix& A, Matrix& At, const double tol, const uint16_t maxit, const bool verbose, const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w) {

    // FIXME: segfaults 
    Eigen::MatrixXd h(w.rows(), A.cols());
    Eigen::VectorXd d(w.rows());
    double tol_ = 1;
    if (verbose)
        Rprintf("\n%4s | %8s \n---------------\n", "iter", "tol");

    // alternating least squares update loop
    for (uint16_t iter_ = 0; iter_ < maxit && tol_ > tol; ++iter_) {

        Eigen::MatrixXd w_it = w;
        // update h
        predict(A, w, h, L1, L2, threads);
        scale(h, d);
        Rcpp::checkUserInterrupt();
        // update w
        predict(At, h, w, L1, L2, threads);
        scale(w, d);

        // calculate tolerance of the model fit to detect convergence
        // correlation between "w" across consecutive iterations
        tol_ = cor(w, w_it);

        if (verbose)
            Rprintf("%4d | %8.2e\n", iter_ + 1, tol_);
        Rcpp::checkUserInterrupt();

    }

    return Rcpp::List::create(Rcpp::Named("w") = w, Rcpp::Named("d") = d, Rcpp::Named("h") = h);

}

//[[Rcpp::export]]
Rcpp::List c_nmf(Rcpp::SparseMatrix& A, Rcpp::SparseMatrix& At, const double tol, const uint16_t maxit, const bool verbose,
                 const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w) {
    return c_nmf_base(A, At, tol, maxit, verbose, L1, L2, threads, w);
}

//[[Rcpp::export]]
Rcpp::List c_nmf_dense(Eigen::MatrixXd& A, Eigen::MatrixXd& At, const double tol, const uint16_t maxit, const bool verbose, const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w) {
    return c_nmf_base(A, At, tol, maxit, verbose, L1, L2, threads, w);
}

// NMF FUNCTION ---------------------------------------------------------------------------------------
// run NMF with linking
//[[Rcpp::export]]
Rcpp::List c_linked_nmf(Rcpp::SparseMatrix A, Rcpp::SparseMatrix At, const double tol, const uint16_t maxit, const bool verbose,
                        const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w, Eigen::MatrixXd link_h,
                        Eigen::MatrixXd link_w) {
    Eigen::MatrixXd h(w.rows(), A.cols());
    Eigen::VectorXd d(w.rows());
    const bool linking_h = (link_h.cols() == A.cols());
    const bool linking_w = (link_w.cols() == A.rows());
    double tol_ = 1;

    if (verbose)
        Rprintf("\n%4s | %8s \n---------------\n", "iter", "tol");

    // alternating least squares update loop

    for (uint16_t iter_ = 0; iter_ < maxit && tol_ > tol; ++iter_) {
        Eigen::MatrixXd w_it = w;
        linking_h ? predict_link(A, w, h, L1, L2, threads, link_h) : predict(A, w, h, L1, L2, threads);
        scale(h, d);
        Rcpp::checkUserInterrupt();
        linking_w ? predict_link(At, h, w, L1, L2, threads, link_w) : predict(At, h, w, L1, L2, threads);
        scale(w, d);
        tol_ = cor(w, w_it);
        if (verbose)
            Rprintf("%4d | %8.2e\n", iter_ + 1, tol_);
        Rcpp::checkUserInterrupt();
    }
    return Rcpp::List::create(Rcpp::Named("w") = w, Rcpp::Named("d") = d, Rcpp::Named("h") = h);
}

// NMF FUNCTION ---------------------------------------------------------------------------------------
// run NMF with masking of a random speckled test set
template <class Matrix>
Rcpp::List c_ard_nmf_base(Matrix& A, Matrix& At, const double tol, const uint16_t maxit, const bool verbose,
                          const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w, const uint32_t seed,
                          const uint32_t inv_density, const double overfit_threshold, const uint16_t trace_test_mse) {
    Eigen::MatrixXd h(w.rows(), A.cols());
    Eigen::VectorXd d(w.rows());
    double tol_ = 1;
    Rcpp::NumericVector test_mse, fit_tol, score_overfit;
    Rcpp::IntegerVector iter;

    if (verbose)
        Rprintf("\n%4s | %8s | %8s \n---------------------------\n", "iter", "tol", "overfit");

    // alternating least squares update loop
    uint16_t iter_ = 0;
    for (; iter_ < maxit && tol_ > tol; ++iter_) {
        Eigen::MatrixXd w_it = w;
        predict_mask(A, seed, inv_density, w, h, L1, L2, threads, false);
        scale(h, d);
        Rcpp::checkUserInterrupt();
        predict_mask(At, seed, inv_density, h, w, L1, L2, threads, true);
        scale(w, d);
        tol_ = cor(w, w_it);
        if (verbose) Rprintf("%4d | %8.2e", iter_ + 1, tol_);
        if (iter_ % trace_test_mse == 0) {
            test_mse.push_back(mse_test(A, w, d, h, seed, inv_density, threads));
            iter.push_back(iter_);
            fit_tol.push_back(tol_);
            double this_err = test_mse(test_mse.size() - 1);
            double min_err = Rcpp::min(test_mse);
            score_overfit.push_back((this_err - min_err) / (this_err + min_err));
            if (verbose) Rprintf(" | %8.2e\n", score_overfit[score_overfit.size() - 1]);
            if (score_overfit[score_overfit.size() - 1] > overfit_threshold)
                break;
        } else if (verbose)
            Rprintf(" | %8s\n", "-");
        Rcpp::checkUserInterrupt();
    }
    if (iter_ % trace_test_mse != 0) {
        test_mse.push_back(mse_test(A, w, d, h, seed, inv_density, threads));
        iter.push_back(iter_);
        fit_tol.push_back(tol_);
        if (test_mse.size() > 0) {
            double min_err = Rcpp::min(test_mse);
            double this_err = test_mse(test_mse.size() - 1);
            score_overfit.push_back((this_err - min_err) / (this_err + min_err));
        } else {
            score_overfit.push_back(0.0);
        }
    }

    // calculate test set reconstruction error, if applicable
    return Rcpp::List::create(
        Rcpp::Named("w") = w,
        Rcpp::Named("d") = d,
        Rcpp::Named("h") = h,
        Rcpp::Named("test_mse") = test_mse,
        Rcpp::Named("iter") = iter,
        Rcpp::Named("tol") = fit_tol,
        Rcpp::Named("score_overfit") = score_overfit);
}

//[[Rcpp::export]]
Rcpp::List c_ard_nmf(Rcpp::SparseMatrix& A, Rcpp::SparseMatrix& At, const double tol, const uint16_t maxit, const bool verbose,
                     const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w, const uint32_t seed,
                     const uint32_t inv_density, const double overfit_threshold, const uint16_t trace_test_mse) {
    return c_ard_nmf_base(A, At, tol, maxit, verbose, L1, L2, threads, w, seed, inv_density, overfit_threshold, trace_test_mse);
}

//[[Rcpp::export]]
Rcpp::List c_ard_nmf_dense(Eigen::MatrixXd& A, Eigen::MatrixXd& At, const double tol, const uint16_t maxit, const bool verbose,
                           const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w, const uint32_t seed,
                           const uint32_t inv_density, const double overfit_threshold, const uint16_t trace_test_mse) {
    return c_ard_nmf_base(A, At, tol, maxit, verbose, L1, L2, threads, w, seed, inv_density, overfit_threshold, trace_test_mse);
}

//[[Rcpp::export]]
Rcpp::S4 log_normalize(Rcpp::SparseMatrix A_, const unsigned int scale_factor, const int threads) {
    Rcpp::SparseMatrix A = A_.clone();

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < A.cols(); ++i) {
        // calculate column sum
        double sum = 0;
        for (Rcpp::SparseMatrix::InnerIterator it(A, i); it; ++it)
            sum += it.value();

        // multiply by scale_factor
        double norm = scale_factor / sum;

        // log-transform
        for (Rcpp::SparseMatrix::InnerIterator it(A, i); it; ++it)
            it.value() = std::log1p(it.value() * norm);
    }
    return A.wrap();
}
