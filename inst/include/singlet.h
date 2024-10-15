#ifndef SINGLET_H
#define SINGLET_H

#include <RcppCommon.h>

// forward declare Rcpp::as<> Exporter
namespace Rcpp {
class SparseMatrix;
namespace traits {
template <>
class Exporter<Rcpp::SparseMatrix>;
}  // namespace traits
}  // namespace Rcpp

//[[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// now pull in IVSparse after loading Eigen
#include <IVSparse.h>

// this class is provided for consistency with Eigen::SparseMatrix, but using
// R objects (i.e. Rcpp::NumericVector, Rcpp::IntegerVector) that comprise Matrix::dgCMatrix in R.
// R objects are pointers to underlying memory-mapped SEXP vectors, and are usable in C++ without any
// affect on performance. Thus, this class achieves zero-copy access to R sparse matrix objects, with equal
// performance for read-only column iteration (`InnerIterator`) like `Eigen::SparseMatrix<double>`.
//
// The class is designed with an `InnerIterator` class that exactly mimics `Eigen::SparseMatrix<T>::InnerIterator`,
// and also contains `.rows()` and `.cols()` member functions. This allows it to substitute for `Eigen::SparseMatrix`
// in all SLAM routines.
namespace Rcpp {
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
#endif