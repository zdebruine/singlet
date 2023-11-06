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

// update h given A and w
inline void predict(std::vector<Rcpp::SparseMatrix> A, const Eigen::MatrixXd& w, Eigen::MatrixXd& h, const double L1, const double L2, const int threads) {
    Eigen::MatrixXd a = AAt(w);
    if (L2 != 0) a.diagonal().array() *= (1 - L2);
    size_t offset = 0;
    for (size_t chunk = 0; chunk < A.size(); ++chunk) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
        for (size_t i = 0; i < A[chunk].cols(); ++i) {
            if (A[chunk].p[i] == A[chunk].p[i + 1]) continue;
            Eigen::VectorXd b = Eigen::VectorXd::Zero(h.rows());
            for (Rcpp::SparseMatrix::InnerIterator it(A[chunk], i); it; ++it)
                b += it.value() * w.col(it.row());
            b.array() -= L1;
            nnls(a, b, h, i + offset);
        }
        offset += A[chunk].cols();
    }
}

//[[Rcpp::export]]
Rcpp::List c_project_model(Rcpp::SparseMatrix A, Eigen::MatrixXd w, const double L1, const double L2, const int threads) {
    if (w.rows() == A.rows()) w = w.transpose();
    Eigen::VectorXd d = Eigen::VectorXd::Ones(w.rows());
    scale(w, d);
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(w.rows(), A.cols());
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

// update h given A and w
inline void predict_mask(std::vector<Rcpp::SparseMatrix>& A, rng seed, const uint64_t inv_density, const Eigen::MatrixXd& w,
                         Eigen::MatrixXd& h, const double L1, const double L2, const int threads, const bool mask_t) {
    Eigen::MatrixXd a = AAt(w);
    if (L2 != 0) a.diagonal().array() *= (1 - L2);
    size_t offset = 0;
    for (size_t chunk = 0; chunk < A.size(); ++chunk) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
        for (size_t i = 0; i < A[chunk].cols(); ++i) {
            if (A[chunk].p[i] == A[chunk].p[i + 1]) continue;
            Eigen::VectorXd b = Eigen::VectorXd::Zero(h.rows());
            Rcpp::SparseMatrix::InnerIterator it(A[chunk], i);
            std::vector<uint64_t> idx;
            idx.reserve(w.cols() / inv_density);
            for (uint64_t j = 0; j < w.cols(); ++j) {
                if (mask_t ? seed.draw(j, i + offset, inv_density) : seed.draw(i + offset, j, inv_density)) {
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
            nnls(a_i, b, h, i + offset);
        }
        offset += A[chunk].cols();
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
inline double mse_test(std::vector<Rcpp::SparseMatrix> A, const Eigen::MatrixXd& w, Eigen::VectorXd& d, Eigen::MatrixXd& h,
                       rng seed, const uint64_t inv_density, const uint16_t threads) {
    // multiply w by d
    Eigen::MatrixXd w_ = w.transpose();
    for (size_t i = 0; i < w.cols(); ++i)
        for (size_t j = 0; j < w.rows(); ++j)
            w_(i, j) *= d(j);

    Eigen::VectorXd losses(h.cols());
    size_t offset = 0;
    for (size_t chunk = 0; chunk < A.size(); ++chunk) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
        for (size_t j = 0; j < A[chunk].cols(); ++j) {
            uint64_t n = 0;
            double s = 0;
            Rcpp::SparseMatrix::InnerIterator it(A[chunk], j);
            for (uint64_t i = 0; i < w_.rows(); ++i) {
                if (seed.draw(j + offset, i, inv_density)) {
                    ++n;
                    if (it && i == it.row()) {
                        s += std::pow(w_.row(i) * h.col(j + offset) - it.value(), 2);
                        ++it;
                    } else {
                        s += std::pow(w_.row(i) * h.col(j + offset), 2);
                    }
                } else if (it && i == it.row()) {
                    ++it;
                }
            }
            losses(j + offset) = (n > 0) ? s / n : 0;
        }
        offset += A[chunk].cols();
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
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(w.rows(), A.cols());
    Eigen::VectorXd d = Eigen::VectorXd::Ones(w.rows());
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

// NMF FUNCTIONS ---------------------------------------------------------------------------------------
// no linking or masking
//[[Rcpp::export]]
Rcpp::List c_nmf_sparse_list(Rcpp::List A_, Rcpp::List& At_, const double tol, const uint16_t maxit, const bool verbose, const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w) {
    std::vector<Rcpp::SparseMatrix> A, At;
    for (auto& A_i : A_)
        A.push_back(A_i);
    for (auto& At_i : At_)
        At.push_back(At_i);

    size_t m = A[0].rows(), n = At[0].rows();

    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(w.rows(), n);
    Eigen::VectorXd d = Eigen::VectorXd::Ones(w.rows());
    double tol_ = 1;
    if (verbose)
        Rprintf("\n%4s | %8s \n---------------\n", "iter", "tol");

    // alternating least squares update loop
    for (uint16_t iter_ = 0; iter_ < maxit && tol_ > tol; ++iter_) {
        Eigen::MatrixXd w_it = w;
        predict(A, w, h, L1, L2, threads);
        scale(h, d);
        Rcpp::checkUserInterrupt();
        predict(At, h, w, L1, L2, threads);
        scale(w, d);
        tol_ = cor(w, w_it);
        if (verbose) Rprintf("%4d | %8.2e\n", iter_ + 1, tol_);
        Rcpp::checkUserInterrupt();
    }
    return Rcpp::List::create(Rcpp::Named("w") = w, Rcpp::Named("d") = d, Rcpp::Named("h") = h);
}

inline void predict(Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& w, Eigen::MatrixXd& h, const double L1, const double L2, const int threads) {
    Eigen::MatrixXd a = AAt(w);
    if (L2 != 0) a.diagonal().array() *= (1 - L2);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < h.cols(); ++i) {
        // TO DO: check for empty column
        Eigen::VectorXd b = Eigen::VectorXd::Zero(h.rows());
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it)
            b += it.value() * w.col(it.row());
        b.array() -= L1;
        nnls(a, b, h, i);
    }
}

// NMF FUNCTIONS ---------------------------------------------------------------------------------------
// no linking or masking

// mega-kudos to the amazing Bing GPT-4 chat engine which got this entire function on the first try
// A function that takes two Eigen::SparseMatrix objects as input and returns a new Eigen::SparseMatrix object that is the column-wise concatenation of the inputs
//[[Rcpp::export]]
Eigen::SparseMatrix<double> cbind_Eigen(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B) {
    // Check that the inputs have the same number of rows
    assert(A.rows() == B.rows());
    // Create a new sparse matrix with the same number of rows and the sum of the number of columns of the inputs
    Eigen::SparseMatrix<double> C(A.rows(), A.cols() + B.cols());
    // Reserve enough space for the non-zero elements
    C.reserve(A.nonZeros() + B.nonZeros());
    // Loop over the columns of the inputs
    for (size_t c = 0; c < A.cols() + B.cols(); ++c) {
        // Start a new column in the output matrix
        C.startVec(c);
        // If the current column belongs to the first input matrix
        if (c < A.cols()) {
            // Copy the non-zero elements from the first input matrix
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, c); it; ++it) {
                C.insertBack(it.row(), c) = it.value();
            }
        }
        // Otherwise, if the current column belongs to the second input matrix
        else {
            // Copy the non-zero elements from the second input matrix
            for (Eigen::SparseMatrix<double>::InnerIterator it(B, c - A.cols()); it; ++it) {
                C.insertBack(it.row(), c) = it.value();
            }
        }
    }
    // Finalize the output matrix
    C.finalize();
    // Return the output matrix
    return C;
}

// A function that takes an Rcpp::List of dgCMatrix as input and returns a list of Eigen::SparseMatrix
//[[Rcpp::export]]
Eigen::SparseMatrix<double> convert_dgCMatrix_to_SparseMatrix(const Rcpp::List& L, const bool verbose = true) {
    // Get the length of the input list
    int n = L.size();
    // Create an output list of the same length
    Rcpp::List out(n);
    // Loop over the elements of the input list
    for (int i = 0; i < n; i++) {
        // Get the current element as an S4 object
        Rcpp::S4 mat = L[i];
        // Get the i, x, and p vectors from the S4 object
        Rcpp::IntegerVector i_vec = mat.slot("i");
        Rcpp::NumericVector x_vec = mat.slot("x");
        Rcpp::IntegerVector p_vec = mat.slot("p");
        // Get the dimensions of the matrix from the S4 object
        Rcpp::IntegerVector dim = mat.slot("Dim");
        int rows = dim[0];
        int cols = dim[1];
        // Create an Eigen::SparseMatrix object with the same dimensions
        Eigen::SparseMatrix<double> sp_mat(rows, cols);
        // Reserve enough space for the non-zero elements
        sp_mat.reserve(x_vec.size());
        // Loop over the columns of the matrix
        for (int j = 0; j < cols; j++) {
            // Start a new column in the sparse matrix
            sp_mat.startVec(j);
            // Get the range of non-zero elements in the current column
            int start = p_vec[j];
            int end = p_vec[j + 1];
            // Loop over the non-zero elements in the current column
            for (int k = start; k < end; k++) {
                // Insert the value and the row index in the sparse matrix
                sp_mat.insertBack(i_vec[k], j) = x_vec[k];
            }
        }
        // Finalize the sparse matrix
        sp_mat.finalize();
        // Store the sparse matrix in the output list
        out[i] = sp_mat;
    }
    if (verbose) Rprintf("preparing to cbind matrices\n");

    // recursively cbind Eigen::SparseMatrix objects in the list
    if (verbose) Rprintf("cbinding matrix 2/%i\n", out.size());
    Eigen::SparseMatrix<double> result = cbind_Eigen(out[0], out[1]);
    if (verbose) Rprintf("   ncol = %i\n", result.cols());
    for (size_t i = 2; i < out.size(); ++i) {
        if (verbose) Rprintf("cbinding matrix %i/%i\n", i + 1, out.size());
        result = cbind_Eigen(result, out[i]);
        if (verbose) Rprintf("   ncol = %i\n", result.cols());
    }
    return result;
}

//[[Rcpp::export]]
Rcpp::List run_nmf_on_dgCMatrix_list(Rcpp::List A_, const double tol, const uint16_t maxit, const bool verbose, const uint16_t threads, Eigen::MatrixXd w) {
    Eigen::SparseMatrix<double> A = convert_dgCMatrix_to_SparseMatrix(A_);
    Eigen::SparseMatrix<double> At = A.transpose();

    if (w.cols() != A.rows()) Rcpp::stop("number of rows in 'w' and 'A' is incompatible!");
    size_t m = A.rows(), n = At.rows();

    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(w.rows(), n);
    Eigen::VectorXd d = Eigen::VectorXd::Ones(w.rows());
    double tol_ = 1;
    if (verbose)
        Rprintf("\n%4s | %8s \n---------------\n", "iter", "tol");

    // alternating least squares update loop
    for (uint16_t iter_ = 0; iter_ < maxit && tol_ > tol; ++iter_) {
        Eigen::MatrixXd w_it = w;
        predict(A, w, h, 0, 0, threads);
        scale(h, d);
        Rcpp::checkUserInterrupt();
        predict(At, h, w, 0, 0, threads);
        scale(w, d);
        tol_ = cor(w, w_it);
        if (verbose) Rprintf("%4d | %8.2e\n", iter_ + 1, tol_);
        Rcpp::checkUserInterrupt();
    }
    return Rcpp::List::create(Rcpp::Named("w") = w, Rcpp::Named("d") = d, Rcpp::Named("h") = h);
}

// NMF FUNCTIONS ---------------------------------------------------------------------------------------
// no linking or masking
template <class Matrix>
Rcpp::List c_mu_nmf_base(Matrix& A, Matrix& At, const double tol, const uint16_t maxit, const bool verbose, const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w) {
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(w.rows(), A.cols());
    h = h.array().abs();
    Eigen::VectorXd d = Eigen::VectorXd::Ones(w.rows());
    double tol_ = 1;
    if (verbose)
        Rprintf("\n%4s | %8s \n---------------\n", "iter", "tol");

    // alternating least squares update loop
    for (uint16_t iter_ = 0; iter_ < maxit && tol_ > tol; ++iter_) {
        Eigen::MatrixXd w_it = w;
        // update h

        Eigen::VectorXd w_rowsums = w.rowwise().squaredNorm();
        for (size_t i = 0; i < h.cols(); ++i) {
            Eigen::VectorXd numer = Eigen::VectorXd::Zero(w.rows());
            Eigen::VectorXd denom = Eigen::VectorXd::Zero(w.rows());
            for (typename Matrix::InnerIterator it(A, i); it; ++it)
                numer += it.value() * w.col(it.row());
            for (size_t j = 0; j < w.rows(); ++j)
                h(j, i) = numer(j) / (w_rowsums(j) * h(j, i));
        }

        Rcpp::checkUserInterrupt();
        // update w
        Eigen::VectorXd h_rowsums = h.rowwise().squaredNorm();
        for (size_t i = 0; i < w.cols(); ++i) {
            Eigen::VectorXd numer = Eigen::VectorXd::Zero(w.rows());
            for (typename Matrix::InnerIterator it(At, i); it; ++it)
                numer += it.value() * h.col(it.row());
            for (size_t j = 0; j < w.rows(); ++j)
                w(j, i) = numer(j) / (h_rowsums(j) * w(j, i));
        }

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
Rcpp::List c_mu_nmf(Rcpp::SparseMatrix& A, Rcpp::SparseMatrix& At, const double tol, const uint16_t maxit, const bool verbose,
                    const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w) {
    return c_mu_nmf_base(A, At, tol, maxit, verbose, L1, L2, threads, w);
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
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(w.rows(), A.cols());
    Eigen::VectorXd d = Eigen::VectorXd::Ones(w.rows());
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
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(w.rows(), A.cols());
    Eigen::VectorXd d = Eigen::VectorXd::Ones(w.rows());
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
Rcpp::List c_ard_nmf_sparse_list(Rcpp::List A_, Rcpp::List At_, const double tol, const uint16_t maxit, const bool verbose,
                                 const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w, const uint64_t rng_seed,
                                 const uint64_t inv_density, const double overfit_threshold, const uint16_t trace_test_mse) {
    std::vector<Rcpp::SparseMatrix> A, At;
    for (auto& A_i : A_)
        A.push_back(A_i);
    for (auto& At_i : At_)
        At.push_back(At_i);

    size_t m = A[0].rows(), n = At[0].rows();

    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(w.rows(), n);
    Eigen::VectorXd d = Eigen::VectorXd::Ones(w.rows());
    h.setZero();
    d.setOnes();
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
Rcpp::List c_ard_nmf_dense(Eigen::MatrixXd& A, Eigen::MatrixXd& At, const double tol, const uint16_t maxit, const bool verbose,
                           const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w, const uint64_t seed,
                           const uint64_t inv_density, const double overfit_threshold, const uint16_t trace_test_mse) {
    return c_ard_nmf_base(A, At, tol, maxit, verbose, L1, L2, threads, w, seed, inv_density, overfit_threshold, trace_test_mse);
}

// ---- CONVOLUTIONAL NMF FUNCTIONS

//[[Rcpp::export]]
Rcpp::S4 spatial_graph(std::vector<double> c1, std::vector<double> c2, double max_dist, size_t max_k = 100, const size_t threads = 0) {
    // calculate euclidean distances between all elements in x and y, return only those less than max_dist
    size_t n = c1.size();
    double scale_factor = 1 / max_dist;
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
        // normalize columns to sum to 1
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

template <typename T>
std::vector<T> apply_permutation(
    const std::vector<T>& vec,
    const std::vector<std::size_t>& p) {
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
                   [&](std::size_t i) { return vec[i]; });
    return sorted_vec;
}

template <typename T>
T jaccard_distance(T* p, T* q, size_t k) {
    T pq = 0, pp = 0, qq = 0;
    for (size_t i = 0; i < k; ++i, ++p, ++q) {
        pq += *p * *q;
        pp += *p * *p;
        qq += *q * *q;
    }
    return 1 - pq / (pp + qq - pq);
}

template <typename T>
T euclidean_distance(T* p, T* q, size_t k) {
    T pq = 0;
    for (size_t i = 0; i < k; ++i, ++p, ++q)
        pq += (*p - *q) * (*p - *q);
    return std::sqrt(pq);
}

template <typename T>
T manhattan_distance(T* p, T* q, size_t k) {
    T pq = 0;
    for (size_t i = 0; i < k; ++i, ++p, ++q)
        pq += std::abs(*p - *q);
    return std::sqrt(pq);
}

template <typename T>
T hamming_distance(T* p, T* q, size_t k) {
    T sum = 0;
    for (size_t i = 0; i < k; ++i, ++p, ++q)
        if (*p != *q) ++sum;
    return sum;
}

template <typename T>
T kullback_distance(T* p, T* q, size_t k) {
    T pdivq = 0, psum = 0;
    for (size_t i = 0; i < k; ++i, ++p, ++q) {
        if (*q != 0) pdivq += *p / *q;
        psum += *p;
    }
    return psum * std::log(pdivq);
}

template <typename T>
T cosine_distance(T* p, T* q, size_t k) {
    T pq = 0, pp = 0, qq = 0;
    for (size_t i = 0; i < k; ++i, ++p, ++q) {
        pq += *p * *q;
        pp += *p * *p;
        qq += *q * *q;
    }
    return 1 - (pq / (std::sqrt(pp) * std::sqrt(qq)));
}

// develop benchmarking objectives that discriminate histology metadata
// Calculate LKNN graph (k, radius, metric, max_dist) then run CNMF at rank of normal NMF
// - develop method for edge detection (cells next to one another that are different in the same ways), all pairwise relative differences in gene expression, run NMF, G, cNMF.
// - measure the contiguousness of a pattern in one dataset, and see how contiguous it is when projected onto another dataset

// Local K-Nearest Neighbors
// Find the K-nearest neighbors on matrix "m" within some radius of one another (on coordinates "coord_x", "coord_y")
// Also impose a minimum distance cutoff to remove very poor-quality nearest neighbors. This cutoff is not particularly useful for Euclidean distance.
//[[Rcpp::export]]
Rcpp::S4 c_LKNN(Eigen::MatrixXf m, Eigen::VectorXf coord_x, Eigen::VectorXf coord_y, size_t k, float radius, std::string metric, bool similarity, float max_dist, bool verbose, size_t threads) {
    if (m.cols() != m.rows() && m.rows() == coord_x.size()) m = m.transpose();
    if (m.cols() != coord_x.size()) Rcpp::stop("number of columns in 'm' must be equal to number of coordinates");
    if (coord_x.size() != coord_y.size()) Rcpp::stop("length of coordinate vectors must be equivalent");

    size_t n_max_edges = std::ceil(std::pow(radius * 2 + 1, 2)) - 1;
    size_t n = m.cols();
    if (verbose) Rprintf("number of edges per node: %u\n", n_max_edges);

    // allocate sufficient memory in the sparse matrix structure
    Eigen::VectorXf x = Eigen::VectorXf::Zero(n_max_edges * n);
    Eigen::VectorXi i = Eigen::VectorXi::Zero(n_max_edges * n);
    Eigen::VectorXi p = Eigen::VectorXi::Zero(n + 1);

    // column pointers assume that the maximum number of edges are used
    for (size_t i = 1; i < n + 1; ++i)
        p(i) = p(i - 1) + n_max_edges;

    if (verbose) Rprintf("filtering %llu edges\n", n * n_max_edges);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t point1 = 0; point1 < n; ++point1) {
        // calculate Jaccard overlap between h[, col] and all other columns within radius
        // TO DO:  Add other distance measures
        std::vector<size_t> i_point1;
        std::vector<float> x_point1;
        for (size_t point2 = 0; point2 < n; ++point2) {
            if (point1 == point2) continue;
            // check if point2 is within radius of point1
            float d = std::sqrt((coord_x[point1] - coord_x[point2]) * (coord_x[point1] - coord_x[point2]) + (coord_y[point1] - coord_y[point2]) * (coord_y[point1] - coord_y[point2]));
            if (d <= radius) {
                // calculate distance between both points in "m"
                float d12;
                if (metric == "jaccard") {
                    d12 = jaccard_distance(&m(0, point1), &m(0, point2), m.rows());
                    if (!similarity) d12 = 1 - d12;
                } else if (metric == "cosine") {
                    d12 = cosine_distance(&m(0, point1), &m(0, point2), m.rows());
                    if (!similarity) d12 = 1 - d12;
                } else if (metric == "manhattan") {
                    d12 = manhattan_distance(&m(0, point1), &m(0, point2), m.rows());
                } else if (metric == "hamming") {
                    d12 = hamming_distance(&m(0, point1), &m(0, point2), m.rows());
                } else if (metric == "kl") {
                    d12 = kullback_distance(&m(0, point1), &m(0, point2), m.rows());
                } else {
                    d12 = euclidean_distance(&m(0, point1), &m(0, point2), m.rows());
                }
                if (max_dist != 0 && d12 > max_dist) continue;
                x_point1.push_back(d12);
                i_point1.push_back(point2);
            }
        }
        if (x_point1.size() > k) {
            // sort i_point1 and x_point1 by value in x, increasing
            std::vector<size_t> sort_perm(x_point1.size(), 0);
            for (size_t j = 0; j < sort_perm.size(); ++j) sort_perm[j] = j;
            std::sort(sort_perm.begin(), sort_perm.end(), [&](const size_t& a, const size_t& b) { return (x_point1[a] < x_point1[b]); });
            x_point1 = apply_permutation(x_point1, sort_perm);
            i_point1 = apply_permutation(i_point1, sort_perm);

            // resize x_point1 and i_point1 to k + 1
            x_point1.resize(k);
            i_point1.resize(k);

            // sort i_point1 and x_point1 by value in i, increasing
            sort_perm.resize(k);
            for (size_t j = 0; j < k; ++j) sort_perm[j] = j;
            std::sort(sort_perm.begin(), sort_perm.end(), [&](const size_t& a, const size_t& b) { return (i_point1[a] < i_point1[b]); });
            x_point1 = apply_permutation(x_point1, sort_perm);
            i_point1 = apply_permutation(i_point1, sort_perm);
        }
        // write x_point1 and i_point1 to x and i vectors
        size_t pos = p[point1];
        for (size_t pos = p[point1], pos2 = 0; pos2 < x_point1.size(); ++pos, ++pos2) {
            x(pos) = x_point1[pos2];
            i(pos) = i_point1[pos2];
        }
    }
    // drop0's
    size_t pos1 = 0, pos2 = 0;
    for (size_t pos_p = 1; pos2 < x.size(); ++pos2) {
        if (pos2 == p(pos_p)) {
            p(pos_p) = pos1;
            ++pos_p;
        }
        if (x(pos2) != 0) {
            x(pos1) = x(pos2);
            i(pos1) = i(pos2);
            ++pos1;
        }
    }
    p(p.size() - 1) = pos1;
    // x.resize(pos1);
    // i.resize(pos1);

    if (verbose) Rprintf("selected %llu edges\n", pos1);

    // wrap Eigen vectors to dgCMatrix
    Rcpp::IntegerVector Dim(2, n);
    Rcpp::IntegerVector i_(pos1), p_(p.size());
    Rcpp::NumericVector x_(pos1);
    for (size_t j = 0; j < pos1; ++j) {
        i_(j) = i(j);
        x_(j) = x(j);
    }
    for (size_t j = 0; j < p.size(); ++j)
        p_(j) = p(j);
    Rcpp::SparseMatrix res(x_, i_, p_, Dim);
    return res.wrap();
}

//[[Rcpp::export]]
Rcpp::S4 c_SNN(Rcpp::SparseMatrix G, double min_similarity, size_t threads) {
    size_t n = G.cols();
    std::vector<size_t> idx, ptrs(n + 1), nnz(n);
    std::vector<double> vals;
    // compute number of nonzeros in each column
    for (size_t i = 0; i < n; ++i) nnz[i] = G.p[i + 1] - G.p[i];

    // TO DO:  parallelize with OpenMP
    // had some challenges getting memory management set up correctly
    // seems like Massimiliano's answer here may have some merit:
    //  https://stackoverflow.com/questions/24765180/parallelizing-a-for-loop-using-openmp-replacing-push-back
    for (size_t i = 0; i < n; ++i) {
        if (nnz[i] != 0) {
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    idx.push_back(i);
                    vals.push_back(1);
                } else {
                    if (nnz[j] != 0) {
                        int* ptr1 = &G.i[G.p[i]];
                        int* ptr1_end = &G.i[G.p[i + 1] - 1];
                        int* ptr2 = &G.i[G.p[j]];
                        int* ptr2_end = &G.i[G.p[j + 1] - 1];
                        size_t intersection = 0;
                        while (ptr1 <= ptr1_end && ptr2 <= ptr2_end) {
                            if (*ptr1 == *ptr2) {
                                ++intersection;
                                ++ptr1;
                                ++ptr2;
                            } else if (*ptr1 < *ptr2)
                                ++ptr1;
                            else
                                ++ptr2;
                        }
                        if (intersection != 0) {
                            double sim = (double)intersection / (double)(nnz[i] + nnz[j] - intersection);
                            if (sim > min_similarity) {
                                idx.push_back(j);
                                vals.push_back(sim);
                            }
                        }
                    }
                }
            }
        }
        ptrs[i + 1] = idx.size();
    }
    Rcpp::NumericVector x_(vals.size());
    Rcpp::IntegerVector i_(idx.size());
    Rcpp::IntegerVector p_(ptrs.size());
    for (size_t j = 0; j < vals.size(); ++j) {
        x_(j) = vals[j];
        i_(j) = idx[j];
    }
    for (size_t j = 0; j < n + 1; ++j)
        p_(j) = ptrs[j];
    Rcpp::IntegerVector Dim(2, G.cols());
    Rcpp::SparseMatrix res(x_, i_, p_, Dim);
    return res.wrap();
}

// update h given A and w
inline void gcnmf_update_h(Rcpp::SparseMatrix A, const Eigen::MatrixXd& w, Eigen::MatrixXd& h, Rcpp::SparseMatrix G, const double L1, const double L2, const int threads) {
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
inline void gcnmf_update_w(Rcpp::SparseMatrix At, Eigen::MatrixXd& w, const Eigen::MatrixXd& h, Rcpp::SparseMatrix G, const double L1, const double L2, const int threads) {
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

//[[Rcpp::export]]
Rcpp::List c_gcnmf(Rcpp::SparseMatrix& A, Rcpp::SparseMatrix& At, Rcpp::SparseMatrix& G, const double tol, const uint16_t maxit, const bool verbose, const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w) {
    if (w.rows() == A.rows() && w.rows() != w.cols()) w = w.transpose();
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(w.rows(), A.cols());
    Eigen::VectorXd d = Eigen::VectorXd::Ones(w.rows());
    double tol_ = 1;
    if (verbose) Rprintf("\n%4s | %8s \n---------------\n", "iter", "tol");
    for (uint16_t iter_ = 0; iter_ < maxit && tol_ > tol; ++iter_) {
        Eigen::MatrixXd w_it = w;
        gcnmf_update_h(A, w, h, G, L1, L2, threads);
        scale(h, d);
        gcnmf_update_w(At, w, h, G, L1, L2, threads);
        scale(w, d);
        Rcpp::checkUserInterrupt();
        tol_ = cor(w, w_it);
        if (verbose) Rprintf("%4d | %8.2e\n", iter_ + 1, tol_);
    }
    return Rcpp::List::create(Rcpp::Named("w") = w.transpose(), Rcpp::Named("d") = d, Rcpp::Named("h") = h);
}

//[[Rcpp::export]]
Rcpp::NumericMatrix c_differentiate_model(Rcpp::NumericMatrix& h, Rcpp::SparseMatrix& G) {
    // for all nodes in G, calculate change in h
    if (h.rows() == G.cols() && h.rows() != h.cols()) h = Rcpp::transpose(h);
    if (h.cols() != G.cols()) Rcpp::stop("dimensions of 'h' and 'G' are not compatible");
    if (G.i[0] == 0) Rcpp::stop("this graph should not have on-diagonal ones");
    size_t n = h.rows();
    Rcpp::NumericMatrix h_diff(h.rows() * 2, G.cols());
    for (size_t col1 = 0, pos = 0; col1 < G.cols(); ++col1) {
        for (Rcpp::SparseMatrix::InnerIterator it(G, col1); it; ++it, ++pos) {
            for (size_t k = 0, k2 = n; k < n; ++k, ++k2) {
                double diff = h(k, col1) - h(k, it.row());
                if (diff > 0)
                    h_diff(k, pos) = diff;
                else
                    h_diff(k2, pos) = -diff;
            }
        }
    }
    return h_diff;
}

//[[Rcpp::export]]
Rcpp::IntegerMatrix c_assign_cells_to_edge_clusters(Rcpp::SparseMatrix G, Rcpp::IntegerVector h_diff_clusters) {
    // get number of clusters
    size_t num_clusters = 0;
    for (auto& cl : h_diff_clusters)
        if (cl > num_clusters) num_clusters = cl;

    // create matrix of unique clusters in rows and samples in columns
    Rcpp::IntegerMatrix res(num_clusters, G.cols());
    for (size_t col1 = 0, pos = 0; col1 < G.cols(); ++col1) {
        for (Rcpp::SparseMatrix::InnerIterator it(G, col1); it; ++it, ++pos) {
            res(h_diff_clusters(pos), col1) += 1;
        }
    }
    return res;
}
