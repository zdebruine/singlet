#define EIGEN_NO_DEBUG
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//[[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

// SPARSE MATRIX CLASS  ---------------------------------------------------------------------------------
class spmat {
public:
  // public member objects
  Rcpp::NumericVector x;
  Rcpp::IntegerVector i, p, Dim;
  
  // constructors
  spmat(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, Rcpp::IntegerVector Dim) : x(x), i(i), p(p), Dim(Dim) {}
  spmat(const Rcpp::S4& s) {
    x = s.slot("x");
    i = s.slot("i");
    p = s.slot("p");
    Dim = s.slot("Dim");
  }
  spmat() {}
  
  size_t rows() { return Dim[0]; }
  size_t cols() { return Dim[1]; }
  
  // const column iterator
  class InnerIterator {
  public:
    InnerIterator(spmat& ptr, int col) : ptr(ptr), col_(col), index(ptr.p[col]), max_index(ptr.p[col + 1]) {}
    operator bool() const { return (index < max_index); }
    InnerIterator& operator++() {
      ++index;
      return *this;
    }
    double& value() { return ptr.x[index]; }
    int row() const { return ptr.i[index]; }
    
  private:
    spmat& ptr;
    int col_, index, max_index;
  };
  
  spmat clone() {
    Rcpp::NumericVector x_ = Rcpp::clone(x);
    Rcpp::IntegerVector i_ = Rcpp::clone(i);
    Rcpp::IntegerVector p_ = Rcpp::clone(p);
    Rcpp::IntegerVector Dim_ = Rcpp::clone(Dim);
    return spmat(x_, i_, p_, Dim_);
  }
  
  Rcpp::S4 wrap() {
    Rcpp::S4 s(std::string("dgCMatrix"));
    s.slot("x") = x;
    s.slot("i") = i;
    s.slot("p") = p;
    s.slot("Dim") = Dim;
    return s;
  }
};

// RANDOM NUMBER GENERATOR  ---------------------------------------------------------------------------------
// two-dimensional linear congruential random number generator using
//   Marsaglia's xorshift32 algorithm together with xoroshiro128++ and a
//   random seed. The generated RNG is non-correlated with sequences in `i` and `j`.
//   The returned logical is a probability of 1/max and is transpose-identical (independent of order of `i` and `j`).
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
  spmat A(A_);
  A.clone();
  
  // get sum of each group
  Rcpp::NumericVector sums(n_groups);
  for (size_t j = 0; j < split_by.size(); ++j) {
    for (spmat::InnerIterator it(A, j); it; ++it) {
      sums[split_by[j]] += it.value();
    }
  }
  
  // calculate multiplication factor for each group
  for (size_t j = 1; j < sums.size(); ++j)
    sums[j] /= sums[0];
  
  // normalize each column using the group sum relative to the first group sum
  for (size_t i = 0; i < split_by.size(); ++i) {
    if (split_by[i] != 0) {
      for (spmat::InnerIterator it(A, i); it; ++it)
        it.value() /= sums[split_by[i]];
    }
  }
  
  return A.wrap();
}

// NMF HELPER FUNCTIONS ---------------------------------------------------------------------------------
// Pearson correlation between two matrices
template <typename T>
inline T cor(Eigen::Matrix<T, -1, -1>& x, Eigen::Matrix<T, -1, -1>& y) {
  T x_i, y_i, sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
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

// fast symmetric matrix multiplication, A * A.transpose() - double
Eigen::MatrixXd AAt(const Eigen::MatrixXd& A) {
  Eigen::MatrixXd AAt = Eigen::MatrixXd::Zero(A.rows(), A.rows());
  AAt.selfadjointView<Eigen::Lower>().rankUpdate(A);
  AAt.triangularView<Eigen::Upper>() = AAt.transpose();
  AAt.diagonal().array() += 1e-15;
  return AAt;
}

// fast symmetric matrix multiplication, A * A.transpose() - float
Eigen::MatrixXf AAt(const Eigen::MatrixXf& A) {
  Eigen::MatrixXf AAt = Eigen::MatrixXf::Zero(A.rows(), A.rows());
  AAt.selfadjointView<Eigen::Lower>().rankUpdate(A);
  AAt.triangularView<Eigen::Upper>() = AAt.transpose();
  AAt.diagonal().array() += 1e-15;
  return AAt;
}

// subset columns from a matrix (deep copy)
// could not figure out how to do better than a deep copy:
//   https://stackoverflow.com/questions/72100483/matrix-multiplication-of-an-eigen-matrix-for-a-subset-of-columns
template <typename T>
Eigen::Matrix<T, -1, -1> submat(const Eigen::Matrix<T, -1, -1>& x, const std::vector<uint32_t>& col_indices) {
  Eigen::Matrix<T, -1, -1> x_(x.rows(), col_indices.size());
  for (size_t i = 0; i < col_indices.size(); ++i)
    x_.col(i) = x.col(col_indices[i]);
  return x_;
}

// scale rows in w (or h) to sum to 1 and put previous rowsums in d
template <typename T>
void scale(Eigen::Matrix<T, -1, -1>& w, Eigen::Matrix<T, 1, -1>& d) {
  d = w.rowwise().sum();
  d.array() += 1e-15;
  for (size_t i = 0; i < w.rows(); ++i)
    for (size_t j = 0; j < w.cols(); ++j)
      w(i, j) /= d(i);
};

// NNLS SOLVER OF THE FORM ax=b ---------------------------------------------------------------------------------
// optimized and modified from github.com/linxihui/NNLM "c_nnls" function
template <typename T>
inline void nnls(Eigen::Matrix<T, -1, -1>& a, Eigen::Matrix<T, 1, -1>& b, Eigen::Matrix<T, -1, -1>& h, const size_t sample) {
  T tol = 1;
  for (uint8_t it = 0; it < 100 && (tol / b.size()) > 1e-8; ++it) {
    tol = 0;
    for (size_t i = 0; i < h.rows(); ++i) {
      T diff = b(i) / a(i, i);
      if (-diff > h(i, sample)) {
        if (h(i, sample) != 0) {
          b -= a.col(i) * -h(i, sample);
          tol = 1;
          h(i, sample) = 0;
        }
      } else if (diff != 0) {
        h(i, sample) += diff;
        b -= a.col(i) * diff;
        tol += std::abs(diff / (h(i, sample) + 1e-15));
      }
    }
  }
}

// NMF PROJECTION FUNCTIONS  ------------------------------------------------------------------------------

// update h given A and w
template <typename T>
void predict(spmat A, const Eigen::Matrix<T, -1, -1>& w, Eigen::Matrix<T, -1, -1>& h, const T L1, const int threads) {
  Eigen::Matrix<T, -1, -1> a = AAt(w);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
  for (size_t i = 0; i < h.cols(); ++i) {
    if (A.p[i] == A.p[i + 1]) continue;
    Eigen::Matrix<T, 1, -1> b = Eigen::Matrix<T, 1, -1>::Zero(h.rows());
    for (spmat::InnerIterator it(A, i); it; ++it)
      b += it.value() * w.col(it.row());
    b.array() -= L1;
    nnls(a, b, h, i);
  }
}

// update h given A and w while masking values in h
template <typename T>
void predict_link(spmat A, const Eigen::Matrix<T, -1, -1>& w, Eigen::Matrix<T, -1, -1>& h, const T L1, const int threads,
                  const Eigen::Matrix<T, -1, -1>& link_h) {
  Eigen::Matrix<T, -1, -1> a = AAt(w);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
  for (size_t i = 0; i < h.cols(); ++i) {
    if (A.p[i] == A.p[i + 1]) continue;
    Eigen::Matrix<T, 1, -1> b = Eigen::Matrix<T, 1, -1>::Zero(h.rows());
    for (spmat::InnerIterator it(A, i); it; ++it)
      b += it.value() * w.col(it.row());
    b.array() -= L1;
    for (size_t j = 0; j < link_h.rows(); ++j)
      b[j] *= link_h(j, i);
    nnls(a, b, h, i);
  }
}

// update h given A and w while masking a random speckled test set
template <typename T>
void predict_mask(spmat A, const uint32_t seed, const uint32_t inv_density, const Eigen::Matrix<T, -1, -1>& w,
                  Eigen::Matrix<T, -1, -1>& h, const T L1, const int threads, const bool mask_t) {
  Eigen::Matrix<T, -1, -1> a = AAt(w);
  
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
  for (size_t i = 0; i < h.cols(); ++i) {
    if (A.p[i] == A.p[i + 1]) continue;
    Eigen::Matrix<T, 1, -1> b = Eigen::Matrix<T, 1, -1>::Zero(h.rows());
    spmat::InnerIterator it(A, i);
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
    Eigen::Matrix<T, -1, -1> wsub = submat(w, idx);
    Eigen::Matrix<T, -1, -1> asub = AAt(wsub);
    Eigen::Matrix<T, -1, -1> a_i = a - asub;
    nnls(a_i, b, h, i);
  }
}

// CALCULATE ERROR OF TEST SET ----------------------------------------------------------------------------

// calculate mean squared error of the model at test set indices only
template <typename T>
T mse_test(spmat A, const Eigen::Matrix<T, -1, -1>& w, Eigen::Matrix<T, 1, -1>& d, Eigen::Matrix<T, -1, -1>& h,
           const uint32_t seed, const uint32_t inv_density, const uint16_t threads) {
  // multiply w by d
  Eigen::Matrix<T, -1, -1> w_ = w.transpose();
  for (size_t i = 0; i < w.cols(); ++i)
    for (size_t j = 0; j < w.rows(); ++j)
      w_(i, j) *= d(j);
  
  Eigen::Matrix<T, 1, -1> losses(h.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
  for (uint32_t j = 0; j < h.cols(); ++j) {
    uint32_t n = 0;
    T s = 0;
    spmat::InnerIterator it(A, j);
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
    losses(j) = s / n;
  }
  return (T)losses.sum() / h.cols();
}

// NMF FUNCTION ---------------------------------------------------------------------------------------

// run NMF either with masking or linking
template <typename T>
Rcpp::List run_nmf(const Rcpp::S4& A_, const Rcpp::S4& At_, const T tol, const uint16_t maxit, const bool verbose,
                   const T L1, const uint16_t threads, Eigen::Matrix<T, -1, -1> w, Eigen::Matrix<T, -1, -1> link_h,
                   Eigen::Matrix<T, -1, -1> link_w, const uint32_t seed, const uint32_t inv_density) {
  if (verbose) Rprintf("\n%4s | %8s \n---------------\n", "iter", "tol");
  spmat A(A_), At(At_);
  Eigen::Matrix<T, -1, -1> h(w.rows(), A.cols());
  Eigen::Matrix<T, 1, -1> d(w.rows());
  const bool masking = (inv_density != 0);
  const bool linking_h = (link_h.cols() == A.cols());
  const bool linking_w = (link_w.cols() == A.rows());
  T tol_ = 1;
  for (uint16_t iter_ = 0; iter_ < maxit && tol_ > tol; ++iter_) {
    Eigen::Matrix<T, -1, -1> w_it = w;
    
    // update h using appropriate prediction method
    if (linking_h) {
      predict_link(A, w, h, L1, threads, link_h);
    } else if (masking) {
      predict_mask(A, seed, inv_density, w, h, L1, threads, false);
    } else {
      predict(A, w, h, L1, threads);
    }
    scale(h, d);
    Rcpp::checkUserInterrupt();
    
    // update w using appropriate prediction method
    if (linking_w) {
      predict_link(At, h, w, L1, threads, link_w);
    } else if (masking) {
      predict_mask(At, seed, inv_density, h, w, L1, threads, true);
    } else {
      predict(At, h, w, L1, threads);
    }
    scale(w, d);
    
    // calculate tolerance of the model fit to detect convergence
    tol_ = cor(w, w_it);  // correlation between "w" across consecutive iterations
    if (verbose) Rprintf("%4d | %8.2e\n", iter_ + 1, tol_);
    Rcpp::checkUserInterrupt();
  }
  
  // calculate test set reconstruction error, if applicable
  if (masking) {
    T test_mse = mse_test(A, w, d, h, seed, inv_density, threads);
    return Rcpp::List::create(
      Rcpp::Named("w") = w,
      Rcpp::Named("d") = d,
      Rcpp::Named("h") = h,
      Rcpp::Named("test_mse") = test_mse);
  }
  return Rcpp::List::create(
    Rcpp::Named("w") = w,
    Rcpp::Named("d") = d,
    Rcpp::Named("h") = h);
}

// Rcpp wrapper to handle "double" or "float" type specialization
//[[Rcpp::export]]
Rcpp::List c_nmf(const Rcpp::S4& A_, const Rcpp::S4& At_, const double tol_, const uint16_t maxit, const bool verbose,
                 const double L1_, const uint16_t threads, Rcpp::NumericMatrix w_, Rcpp::NumericMatrix link_h_, Rcpp::NumericMatrix link_w_,
                 const uint32_t seed, const uint32_t inv_density, const bool float_precision) {
  Rcpp::List model;
  if (float_precision) {
    float tol = (float)tol_;
    float L1 = (float)L1_;
    Eigen::MatrixXf w = Rcpp::as<Eigen::MatrixXf>(w_);
    Eigen::MatrixXf link_h = Rcpp::as<Eigen::MatrixXf>(link_h_);
    Eigen::MatrixXf link_w = Rcpp::as<Eigen::MatrixXf>(link_w_);
    model = run_nmf(A_, At_, tol, maxit, verbose, L1, threads, w, link_h, link_w, seed, inv_density);
  } else {
    Eigen::MatrixXd w = Rcpp::as<Eigen::MatrixXd>(w_);
    Eigen::MatrixXd link_h = Rcpp::as<Eigen::MatrixXd>(link_h_);
    Eigen::MatrixXd link_w = Rcpp::as<Eigen::MatrixXd>(link_w_);
    model = run_nmf(A_, At_, tol_, maxit, verbose, L1_, threads, w, link_h, link_w, seed, inv_density);
  }
  return model;
}