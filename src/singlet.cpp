//[[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

#include <singlet.h>

// xorshift64 Linear Congruential Generator
class rng {
   private:
    uint64_t state;

   public:
    rng(uint64_t state) : state(state) {}

    void advance_state() {
        state ^= state << 19;
        state ^= state >> 7;
        state ^= state << 36;
    }

    uint64_t operator*() const {
        return state;
    }

    uint64_t rand() {
        uint64_t x = state ^ (state << 38);
        x ^= x >> 13;
        x ^= x << 23;
    }

    uint64_t rand(uint64_t i) {
        // advance i
        i ^= i << 19;
        i ^= i >> 7;
        i ^= i << 36;

        // add i to state
        uint64_t x = state + i;

        // advance state
        x ^= x << 38;
        x ^= x >> 13;
        x ^= x << 23;

        return x;
    }

    uint64_t rand(uint64_t i, uint64_t j) {
        uint64_t x = rand(i);

        // advance j
        j ^= j >> 7;
        j ^= j << 23;
        j ^= j >> 8;

        // add j to state
        x += j;

        // advance state
        x ^= x >> 7;
        x ^= x << 53;
        x ^= x >> 4;

        return x;
    }

    template <typename T>
    T sample(T max_value) {
        return rand() % max_value;
    }

    template <typename T>
    T sample(uint64_t i, T max_value) {
        return rand(i) % max_value;
    }

    template <typename T>
    T sample(uint64_t i, uint64_t j, T max_value) {
        return rand(i, j) % max_value;
    }

    template <typename T>
    bool draw(T probability) {
        return sample(probability) == 0;
    }

    template <typename T>
    bool draw(uint64_t i, T probability) {
        return sample(i, probability) == 0;
    }

    template <typename T>
    bool draw(uint64_t i, uint64_t j, T probability) {
        sample(i, j, probability);
        return sample(i, j, probability) == 0;
    }

    template <typename T>
    double uniform() {
        T x = (T)rand() / UINT64_MAX;
        return x - std::floor(x);
    }

    template <typename T>
    double uniform(uint64_t i) {
        T x = (T)rand(i) / UINT64_MAX;
        return x - std::floor(x);
    }

    template <typename T>
    double uniform(uint64_t i, uint64_t j) {
        T x = (T)rand(i, j) / UINT64_MAX;
        return x - std::floor(x);
    }
};

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
Rcpp::NumericMatrix rowwise_compress_sparse(Rcpp::SparseMatrix& A, const size_t n = 10, const size_t threads = 0) {
    const size_t n_rows = (size_t)std::floor(A.rows() / n);
    Rcpp::NumericMatrix res(n_rows, A.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t col = 0; col < A.cols(); ++col)
        for (Rcpp::SparseMatrix::InnerIterator it(A, col); it; ++it) {
            size_t row = (size_t)std::floor(it.row() / n);
            res(row, col) += it.value();
        }
    for (size_t j = 0; j < res.cols(); ++j)
        for (size_t i = 0; i < res.rows(); ++i)
            res(i, j) /= n;  // calculate mean
    return res;
}

//[[Rcpp::export]]
Rcpp::NumericMatrix rowwise_compress_dense(Rcpp::NumericMatrix& A, const size_t n = 10, const size_t threads = 0) {
    const size_t n_rows = (size_t)std::floor(A.rows() / n);
    Rcpp::NumericMatrix res(n_rows, A.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t col = 0; col < A.cols(); ++col) {
        size_t res_row = 0;
        for (size_t row = 0; row < A.rows(); row += n, ++res_row) {
            for (size_t i = 0; i < n; ++i)
                res(res_row, col) += A(row + i, col);
            res(res_row, col) /= n;  // calculate mean
        }
    }
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
inline Eigen::MatrixXd submat(const Eigen::MatrixXd& x, const std::vector<uint64_t>& col_indices) {
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
        Eigen::VectorXd b = w * A.col(i);
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
inline void predict_mask(Rcpp::SparseMatrix A, rng seed, const uint64_t inv_density, const Eigen::MatrixXd& w,
                         Eigen::MatrixXd& h, const double L1, const double L2, const int threads, const bool mask_t) {
    Eigen::MatrixXd a = AAt(w);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < h.cols(); ++i) {
        if (A.p[i] == A.p[i + 1]) continue;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(h.rows());
        Rcpp::SparseMatrix::InnerIterator it(A, i);
        std::vector<uint64_t> idx;
        idx.reserve(A.rows() / inv_density);
        for (uint64_t j = 0; j < A.rows(); ++j) {
            if (mask_t ? seed.draw(j, i, inv_density) : seed.draw(i, j, inv_density)) {
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
inline void predict_mask(const Eigen::MatrixXd& A, rng seed, const uint64_t inv_density, const Eigen::MatrixXd& w,
                         Eigen::MatrixXd& h, const double L1, const double L2, const int threads, const bool mask_t) {
    Eigen::MatrixXd a = AAt(w);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < h.cols(); ++i) {
        Eigen::VectorXd b = Eigen::VectorXd::Zero(h.rows());
        std::vector<uint64_t> idx;
        idx.reserve(A.rows() / inv_density);
        for (uint64_t j = 0; j < A.rows(); ++j) {
            if (mask_t ? seed.draw(j, i, inv_density) : seed.draw(i, j, inv_density)) {
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
                       rng seed, const uint64_t inv_density, const uint16_t threads) {
    // multiply w by d
    Eigen::MatrixXd w_ = w.transpose();
    for (size_t i = 0; i < w.cols(); ++i)
        for (size_t j = 0; j < w.rows(); ++j)
            w_(i, j) *= d(j);

    Eigen::VectorXd losses(h.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (uint64_t j = 0; j < h.cols(); ++j) {
        uint64_t n = 0;
        double s = 0;
        Rcpp::SparseMatrix::InnerIterator it(A, j);
        for (uint64_t i = 0; i < A.rows(); ++i) {
            if (seed.draw(j, i, inv_density)) {
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
                       rng seed, const uint64_t inv_density, const uint16_t threads) {
    // multiply w by d
    Eigen::MatrixXd w_ = w.transpose();
    for (size_t i = 0; i < w.cols(); ++i)
        for (size_t j = 0; j < w.rows(); ++j)
            w_(i, j) *= d(j);

    Eigen::VectorXd losses(h.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (uint64_t j = 0; j < h.cols(); ++j) {
        uint64_t n = 0;
        double s = 0;
        for (uint64_t i = 0; i < A.rows(); ++i) {
            if (seed.draw(j, i, inv_density)) {
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
                          const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w, const uint64_t rng_seed,
                          const uint64_t inv_density, const double overfit_threshold, const uint16_t trace_test_mse) {
    Eigen::MatrixXd h(w.rows(), A.cols());
    Eigen::VectorXd d(w.rows());
    double tol_ = 1;
    Rcpp::NumericVector test_mse, fit_tol, score_overfit;
    Rcpp::IntegerVector iter;

    rng seed(rng_seed);

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
                     const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w, const uint64_t seed,
                     const uint64_t inv_density, const double overfit_threshold, const uint16_t trace_test_mse) {
    return c_ard_nmf_base(A, At, tol, maxit, verbose, L1, L2, threads, w, seed, inv_density, overfit_threshold, trace_test_mse);
}

//[[Rcpp::export]]
Rcpp::List c_ard_nmf_dense(Eigen::MatrixXd& A, Eigen::MatrixXd& At, const double tol, const uint16_t maxit, const bool verbose,
                           const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w, const uint64_t seed,
                           const uint64_t inv_density, const double overfit_threshold, const uint16_t trace_test_mse) {
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

// ---- CONVOLUTIONAL NMF FUNCTIONS

//[[Rcpp::export]]
Rcpp::S4 spatial_graph(std::vector<double> c1, std::vector<double> c2, double max_dist, size_t max_k = 100, const size_t threads = 0) {
    // calculate euclidean distances between all elements in x and y, return only those less than max_dist
    size_t n = c1.size();
    float scale_factor = 1 / max_dist;
    Eigen::MatrixXd x_ = Eigen::MatrixXd::Zero(max_k, n);
    Eigen::MatrixXi i_ = Eigen::MatrixXi::Zero(max_k, n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < n; ++i) {
        size_t pos = 0;
        for (size_t j = 0; j < n && pos < max_k; ++j) {
            double d = std::sqrt((c1[i] - c1[j]) * (c1[i] - c1[j]) + (c2[i] - c2[j]) * (c2[i] - c2[j]));
            if (d < max_dist) {
                i_(pos, i) = j;
                x_(pos, i) = (max_dist - d) * scale_factor;
                ++pos;
            }
        }
        //normalize columns to sum to 1
        double sum = x_.col(i).sum();
        x_.col(i).array() /= sum;
    }

    // convert to vector form
    size_t nnz = 0;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < max_k; ++j)
            if (x_(j, i) != 0) ++nnz;
    //    for (auto it = i_.data(); it != (i_.data() + i_.size()); ++it)
    //       if (*it != 0) ++nnz;

    Rcpp::IntegerVector idx(nnz), colptrs(n + 1), Dim(2, n);
    Rcpp::NumericVector vals(nnz);
    size_t pos = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < max_k; ++j) {
            if (x_(j, i) != 0) {
                idx[pos] = i_(j, i);
                vals[pos] = x_(j, i);
                ++pos;
            }
        }
        colptrs[i + 1] = pos;
    }
    Rcpp::SparseMatrix m(vals, idx, colptrs, Dim);
    return m.wrap();
}

// update h given A and w
inline void cnmf_update_h(Rcpp::SparseMatrix A, const Eigen::MatrixXd& w, Eigen::MatrixXd& h, Rcpp::SparseMatrix G, const double L1, const double L2, const int threads) {
    Eigen::MatrixXd a = AAt(w);
    if (L2 != 0) a.diagonal().array() *= (1 - L2);
    // calculate b like in non-convolutional updates
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(h.rows(), A.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t j = 0; j < h.cols(); ++j) {
        if (A.p[j] == A.p[j + 1]) continue;
        for (Rcpp::SparseMatrix::InnerIterator it(A, j); it; ++it)
            b.col(j) += it.value() * w.col(it.row());
    }

    // generate convolutional b by convolving b using the graph pattern and then solve NNLS
    for (size_t j = 0; j < b.cols(); ++j) {
        Eigen::VectorXd b_ = Eigen::VectorXd::Zero(b.rows());
        for (Rcpp::SparseMatrix::InnerIterator it(G, j); it; ++it)
            b_ += it.value() * b.col(it.row());
        b_.array() -= L1;
        nnls(a, b_, h, j);
    }
}

// update w given A and h
inline void cnmf_update_w(Rcpp::SparseMatrix At, Eigen::MatrixXd& w, const Eigen::MatrixXd& h, Rcpp::SparseMatrix G, const double L1, const double L2, const int threads) {
    Eigen::MatrixXd a = AAt(h);
    if (L2 != 0) a.diagonal().array() *= (1 - L2);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t j = 0; j < w.cols(); ++j) {
        Eigen::VectorXd b = Eigen::VectorXd::Zero(w.rows());
        for (Rcpp::SparseMatrix::InnerIterator it(At, j); it; ++it) {
            for (Rcpp::SparseMatrix::InnerIterator it2(G, it.row()); it2; ++it2) {
                b += (it.value() * it2.value()) * h.col(it2.row());
            }
        }
        b.array() -= L1;
        nnls(a, b, w, j);
    }
}

// calculate mean squared error of the model at test set indices only
Eigen::VectorXd mse(Rcpp::SparseMatrix A, const Eigen::MatrixXd& w, Eigen::VectorXd& d, Eigen::MatrixXd& h, const uint16_t threads) {
    // multiply w by d
    Eigen::MatrixXd w_ = w.transpose();
    for (size_t i = 0; i < w.cols(); ++i)
        for (size_t j = 0; j < w.rows(); ++j)
            w_(i, j) *= d(j);

    Eigen::VectorXd losses(h.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (uint64_t j = 0; j < h.cols(); ++j) {
        double s = 0;
        Rcpp::SparseMatrix::InnerIterator it(A, j);
        for (uint64_t i = 0; i < A.rows(); ++i) {
            if (it && i == it.row()) {
                s += std::pow(w_.row(i) * h.col(j) - it.value(), 2);
                ++it;
            } else {
                s += std::pow(w_.row(i) * h.col(j), 2);
            }
        }
        losses(j) = s;
    }
    return losses;
}

// update graph
inline void updateGraph(Rcpp::SparseMatrix& G, Rcpp::SparseMatrix& dE) {
    for (size_t i = 0; i < G.cols(); ++i) {
        Rcpp::SparseMatrix::InnerIterator it_dE(dE, i);
        Rcpp::SparseMatrix::InnerIterator it_G(G, i);
        double sum = 0;
        for (; it_G; ++it_G, ++it_dE) {
            it_G.value() = it_G.value() * it_dE.value();
            sum += it_G.value();
        }
        // normalize column to sum to one
        for (Rcpp::SparseMatrix::InnerIterator it(G, i); it; ++it)
            it.value() /= sum;
    }
}

void updateErrGraph(Rcpp::SparseMatrix& dE, Eigen::VectorXd& err, Eigen::VectorXd& prev_err) {
    for (size_t i = 0; i < err.size(); ++i)
        err(i) = (err(i) - prev_err(i)) / prev_err(i);

    for (size_t i = 0; i < dE.cols(); ++i) {
        for (Rcpp::SparseMatrix::InnerIterator it(dE, i); it; ++it)
            it.value() = err[i] - prev_err[i];
    }
}

// we only want pairs of points to be convoluted which have the lowest error.
// If we deconvolute everything, and then convolute everything, find the relative change in Euclidean distance on the "H" matrix of convoluted/deconvoluted and multiply graph weights by that change.

//[[Rcpp::export]]
Rcpp::List c_cnmf(Rcpp::SparseMatrix& A, Rcpp::SparseMatrix& At, Rcpp::SparseMatrix& G, const double tol, const uint16_t maxit, const bool verbose, const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w) {
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(w.rows(), A.cols());
    Eigen::VectorXd d = Eigen::VectorXd::Ones(w.rows());
    double tol_ = 1;
    if (verbose) Rprintf("\n%4s | %8s \n---------------\n", "iter", "tol");
    // store relative change in graph weights from previous iteration (initializing with 1)
    //    Rcpp::SparseMatrix dG = G.clone();
    //    dG.setOnes();
    Rcpp::List costs(maxit);

    // update the graph by least squares to minimize reconstruction error
    //    Eigen::VectorXd err = Eigen::VectorXd::Zero(A.cols());
    //    Eigen::VectorXd prev_err = err;
    //    Eigen::VectorXd rel_err = err;

    // convolve points that the model fits really well to
    // visualize goodness of fit to the model. For a given convolution Gij, the goodness of fit is the euclidean distance between the two points on NMF coordinates.

    for (uint16_t iter_ = 0; iter_ < maxit && tol_ > tol; ++iter_) {
        Eigen::MatrixXd w_it = w;
        cnmf_update_h(A, w, h, G, L1, L2, threads);
        scale(h, d);
        cnmf_update_w(At, w, h, G, L1, L2, threads);
        scale(w, d);
        Rcpp::checkUserInterrupt();

        Rcpp::SparseMatrix cost = G.clone();
        cost.setOnes();
        // calculate Euclidean distance between sample pairs on H matrix
        for (size_t i = 0; i < G.cols(); ++i) {
            for (Rcpp::SparseMatrix::InnerIterator it(cost, i); it; ++it) {
                double d = 0;
                for (size_t k = 0; k < h.rows(); ++k)
                    d += std::pow(h(k, i) - h(k, it.row()), 2);
                it.value() = std::sqrt(d);
            }
        }
        costs[iter_] = cost.wrap();
        // update graph
        //        err = mse(A, w, d, h, threads);  // calculate losses for each column
        //        if (iter_ > 0) {
        //            updateErrGraph(dE, err, prev_err);
        //            updateGraph(G, dE);
        //        }
        //        prev_err = err;

        tol_ = cor(w, w_it);
        if (verbose) Rprintf("%4d | %8.2e\n", iter_ + 1, tol_);
    }
    return Rcpp::List::create(Rcpp::Named("w") = w, Rcpp::Named("d") = d, Rcpp::Named("h") = h, Rcpp::Named("costs") = costs);
}
