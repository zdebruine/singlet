// This file is adapted from RcppEigen headers
//
// It supports R bindings for
//   * Eigen::Matrix<T, -1, -1>
//   * Eigen::Matrix<T, -1, 1>
//
// Copyright (C)      2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
//
// also included is a zero-copy Rcpp SparseMatrix class by Zach DeBruine

#include <RcppCommon.h>

#include "EigenCore"

/* forward declarations */
namespace Rcpp {
namespace traits {
template <typename T>
class Exporter<Eigen::Matrix<T, -1, -1> >;
template <typename T>
class Exporter<Eigen::Matrix<T, -1, 1> >;
}  // namespace traits
}  // namespace Rcpp

namespace Rcpp {
class SparseMatrix;
}  // namespace Rcpp

// forward declare Rcpp::as<> Exporter
namespace Rcpp {
namespace traits {
template <>
class Exporter<Rcpp::SparseMatrix>;
}  // namespace traits
}  // namespace Rcpp

#include <Rcpp.h>

namespace Rcpp {
namespace RcppEigen {

// helper trait to identify if T is a plain object type
template <typename T>
struct is_plain : Rcpp::traits::same_type<T, typename T::PlainObject> {};

// helper trait to identify if the object has dense storage
template <typename T>
struct is_dense : Rcpp::traits::same_type<typename T::StorageKind, Eigen::Dense> {};

// for plain dense objects
template <typename T>
SEXP eigen_wrap_plain_dense(const T& obj, Rcpp::traits::true_type) {
    typename Eigen::internal::conditional<T::IsRowMajor,
                                          Eigen::Matrix<typename T::Scalar,
                                                        T::RowsAtCompileTime,
                                                        T::ColsAtCompileTime>,
                                          const T&>::type objCopy(obj);
    int m = obj.rows(), n = obj.cols();
    R_xlen_t size = static_cast<R_xlen_t>(m) * n;
    SEXP ans = PROTECT(::Rcpp::wrap(objCopy.data(), objCopy.data() + size));
    if (T::ColsAtCompileTime != 1) {
        SEXP dd = PROTECT(::Rf_allocVector(INTSXP, 2));
        int* d = INTEGER(dd);
        d[0] = m;
        d[1] = n;
        ::Rf_setAttrib(ans, R_DimSymbol, dd);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return ans;
}

// plain object, so we can assume data() and size()
template <typename T>
inline SEXP eigen_wrap_is_plain(const T& obj, ::Rcpp::traits::true_type) {
    return eigen_wrap_plain_dense(obj, typename is_dense<T>::type());
}

// when the object is not plain, we need to eval()uate it
template <typename T>
inline SEXP eigen_wrap_is_plain(const T& obj, ::Rcpp::traits::false_type) {
    return eigen_wrap_is_plain(obj.eval(), Rcpp::traits::true_type());
}

template <typename T>
inline SEXP eigen_wrap(const T& obj) {
    return eigen_wrap_is_plain(obj,
                               typename is_plain<T>::type());
}

}  // namespace RcppEigen

// this class is provided for consistency with Eigen::SparseMatrix, but using
// R objects (i.e. Rcpp::NumericVector, Rcpp::IntegerVector) that comprise Matrix::dgCMatrix in R.
// R objects are pointers to underlying memory-mapped SEXP vectors, and are usable in C++ without any
// affect on performance. Thus, this class achieves zero-copy access to R sparse matrix objects, with equal
// performance for read-only column iteration (`InnerIterator`) like `Eigen::SparseMatrix<double>`.
//
// The class is designed with an `InnerIterator` class that exactly mimics `Eigen::SparseMatrix<T>::InnerIterator`,
// and also contains `.rows()` and `.cols()` member functions. This allows it to substitute for `Eigen::SparseMatrix`
// in all SLAM routines.
class SparseMatrix {
public:
  NumericVector x;
  IntegerVector i, p, Dim;
  
  // constructors
  SparseMatrix(NumericVector x, IntegerVector i, IntegerVector p, IntegerVector Dim) : x(x), i(i), p(p), Dim(Dim) {}
  SparseMatrix(const S4& s) {
    if (!s.hasSlot("x") || !s.hasSlot("p") || !s.hasSlot("i") || !s.hasSlot("Dim"))
      throw std::invalid_argument("Cannot construct SparseMatrix from this S4 object");
    x = s.slot("x");
    i = s.slot("i");
    p = s.slot("p");
    Dim = s.slot("Dim");
  }
  SparseMatrix() {}
  
  unsigned int rows() { return Dim[0]; }
  unsigned int cols() { return Dim[1]; }
  
  // const column iterator
  class InnerIterator {
  public:
    InnerIterator(SparseMatrix& ptr, int col) : ptr(ptr), col_(col), index(ptr.p[col]), max_index(ptr.p[col + 1]) {}
    operator bool() const { return (index < max_index); }
    InnerIterator& operator++() {
      ++index;
      return *this;
    }
    double& value() const { return ptr.x[index]; }
    int row() const { return ptr.i[index]; }
    int col() const { return col_; }
    
  private:
    SparseMatrix& ptr;
    int col_, index, max_index;
  };
  
  SparseMatrix clone() {
    NumericVector x_ = Rcpp::clone(x);
    IntegerVector i_ = Rcpp::clone(i);
    IntegerVector p_ = Rcpp::clone(p);
    IntegerVector Dim_ = Rcpp::clone(Dim);
    return SparseMatrix(x_, i_, p_, Dim_);
  }
  
  SparseMatrix transpose() {
    S4 s(std::string("dgCMatrix"));
    s.slot("i") = i;
    s.slot("p") = p;
    s.slot("x") = x;
    s.slot("Dim") = Dim;
    Environment base = Environment::namespace_env("Matrix");
    Function t_r = base["t"];
    S4 At = t_r(_["x"] = s);
    return SparseMatrix(At);
  };
  
  S4 wrap() {
    S4 s(std::string("dgCMatrix"));
    s.slot("x") = x;
    s.slot("i") = i;
    s.slot("p") = p;
    s.slot("Dim") = Dim;
    return s;
  }
};

namespace traits {
/* support for Rcpp::as */

template <typename T, typename value_type>
class MatrixExporterForEigen {
   public:
    typedef value_type r_export_type;

    MatrixExporterForEigen(SEXP x) : object(x) {}
    ~MatrixExporterForEigen() {}

    T get() {
        Shield<SEXP> dims(::Rf_getAttrib(object, R_DimSymbol));
        if (Rf_isNull(dims) || ::Rf_length(dims) != 2) {
            throw ::Rcpp::not_a_matrix();
        }
        int* dims_ = INTEGER(dims);
        T result(dims_[0], dims_[1]);
        value_type* data = result.data();
        ::Rcpp::internal::export_indexing<value_type*, value_type>(object, data);
        return result;
    }

   private:
    SEXP object;
};

// export a dense vector
template <typename T>
class Exporter<Eigen::Matrix<T, -1, 1> >
    : public IndexingExporter<Eigen::Matrix<T, -1, 1>, T> {
   public:
    Exporter(SEXP x) : IndexingExporter<Eigen::Matrix<T, -1, 1>, T>(x) {}
};

// export a dense matrix
template <typename T>
class Exporter<Eigen::Matrix<T, -1, -1> >
    : public MatrixExporterForEigen<Eigen::Matrix<T, -1, -1>, T> {
   public:
    Exporter(SEXP x) : MatrixExporterForEigen<Eigen::Matrix<T, -1, -1>, T>(x) {}
};

// export a sparse matrix
template <>
class Exporter<Rcpp::SparseMatrix> {
  Rcpp::NumericVector x_;
  Rcpp::IntegerVector i, p, Dim;
  
public:
  Exporter(SEXP x) {
    Rcpp::S4 s(x);
    if (!s.hasSlot("x") || !s.hasSlot("p") || !s.hasSlot("i") || !s.hasSlot("Dim"))
      throw std::invalid_argument("Cannot construct Rcpp::SparseMatrix from this S4 object");
    x_ = s.slot("x");
    i = s.slot("i");
    p = s.slot("p");
    Dim = s.slot("Dim");
  }
  
  Rcpp::SparseMatrix get() {
    return Rcpp::SparseMatrix(x_, i, p, Dim);
  }
};

}  // namespace traits
}  // namespace Rcpp