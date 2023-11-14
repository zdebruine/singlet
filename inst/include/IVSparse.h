/**
 * @file SparseMatrix
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief IVSparse Sparse Matrix Library
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

// Library Constants
#define DELIM 0
#define NUM_META_DATA 6
#define META_DATA_SIZE 24
#define ONE_BYTE_MAX 255
#define TWO_BYTE_MAX 65535
#define FOUR_BYTE_MAX 4294967295

// Library Preprocessor Directives

// Parallel Processing Directives (On by default)
#if (defined _OPENMP) && (!defined IVSPARSE_DONT_PARALLEL)
    #define IVSPARSE_HAS_OPENMP
#endif
#ifdef IVSPARSE_HAS_OPENMP
#include <atomic>
#include <omp.h>
#endif

// Debugging Directives (Off by default)
#ifndef IVSPARSE_DEBUG_OFF
#define IVSPARSE_DEBUG
#endif

// Library Includes

// Eigen is already pulled in by "singlet"
//[[Rcpp::depends(RcppEigen)]]
// #include <RcppEigen.h>


#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <type_traits>
#include <iomanip>
#include <type_traits>
// Library Namespaces

// SparseMatrixBase Files
// #include "src/IVSparse_SparseMatrixBase.hpp"
// #include "src/IVSparse_Base_Methods.hpp"

// SparseMatrix Level 3 Files
#include "src/IVCSC/IVCSC_SparseMatrix.hpp"
#include "src/IVCSC/IVCSC_Operators.hpp"
#include "src/IVCSC/IVCSC_Private_Methods.hpp"
#include "src/IVCSC/IVCSC_Methods.hpp"
#include "src/IVCSC/IVCSC_Constructors.hpp"
#include "src/IVCSC/IVCSC_BLAS.hpp"
    // Vector and Iterator Files
    #include "src/Vectors/IVCSC_Vector.hpp"
    #include "src/Vectors/IVCSC_Vector_Methods.hpp"
    #include "src/InnerIterators/IVCSC_Iterator.hpp"
    #include "src/InnerIterators/IVCSC_Iterator_Methods.hpp"

// SparseMatrix Level 2 Files
#include "src/VCSC/VCSC_SparseMatrix.hpp"
#include "src/VCSC/VCSC_Operators.hpp"
#include "src/VCSC/VCSC_Private_Methods.hpp"
#include "src/VCSC/VCSC_Methods.hpp"
#include "src/VCSC/VCSC_Constructors.hpp"
#include "src/VCSC/VCSC_BLAS.hpp"
    // Vector and Iterator Files
    #include "src/Vectors/VCSC_Vector.hpp"
    #include "src/Vectors/VCSC_Vector_Methods.hpp"
    #include "src/InnerIterators/VCSC_Iterator.hpp"
    #include "src/InnerIterators/VCSC_Iterator_Methods.hpp"

// SparseMatrix Level 1 Files
#include "src/CSC/CSC_SparseMatrix.hpp"
#include "src/CSC/CSC_Operators.hpp"
#include "src/CSC/CSC_Private_Methods.hpp"
#include "src/CSC/CSC_Methods.hpp"
#include "src/CSC/CSC_Constructors.hpp"
#include "src/CSC/CSC_BLAS.hpp"
    // Vector and Iterator Files
    #include "src/Vectors/CSC_Vector.hpp"
    #include "src/Vectors/CSC_Vector_Methods.hpp"
    #include "src/InnerIterators/CSC_Iterator.hpp"
    #include "src/InnerIterators/CSC_Iterator_Methods.hpp"