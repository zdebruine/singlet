/**
 * @file CSC_Operators.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Operator Overloads for CSC Sparse Matrices
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

// Assignment Operator
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor> &
SparseMatrix<T, indexT, 1, columnMajor>::operator=(const IVSparse::SparseMatrix<T, indexT, 1, columnMajor> &other) {
  // check for self assignment
  if (this != &other) {
    // free the old data
    if (vals != nullptr) {
      free(vals);
    }
    if (innerIdx != nullptr) {
      free(innerIdx);
    }
    if (outerPtr != nullptr) {
      free(outerPtr);
    }
    if (metadata != nullptr) {
      delete[] metadata;
    }

    // Deep copy the matrix
    metadata = new uint32_t[NUM_META_DATA];
    memcpy(metadata, other.metadata, NUM_META_DATA * sizeof(uint32_t));

    // set the dimensions of the matrix
    numRows = other.numRows;
    numCols = other.numCols;
    outerDim = other.outerDim;
    innerDim = other.innerDim;
    nnz = other.nnz;
    compSize = other.compSize;

    // encode the value type and index type
    encodeValueType();
    index_t = sizeof(indexT);

    // check for an empty matrix
    if (nnz == 0) {
      vals = nullptr;
      innerIdx = nullptr;
      outerPtr = nullptr;
    }

    // allocate the memory
    try {
      vals = (T *)malloc(nnz * sizeof(T));
      innerIdx = (indexT *)malloc(nnz * sizeof(indexT));
      outerPtr = (indexT *)malloc((outerDim + 1) * sizeof(indexT));
    } catch (std::bad_alloc &e) {
      std::cerr << "Error: Failed to allocate memory for the matrix"
                << std::endl;
      exit(1);
    }

    // copy the data
    memcpy(vals, other.vals, nnz * sizeof(T));
    memcpy(innerIdx, other.innerIdx, nnz * sizeof(indexT));
    memcpy(outerPtr, other.outerPtr, (outerDim + 1) * sizeof(indexT));
  }

  // return the matrix
  return *this;
}

// Equality Operator
template <typename T, typename indexT, bool columnMajor>
bool SparseMatrix<T, indexT, 1, columnMajor>::operator==(const SparseMatrix<T, indexT, 1, columnMajor> &other) {
  // check if the dimensions are the same
  if (numRows != other.numRows || numCols != other.numCols) {
    return false;
  }

  // check if the number of nonzeros are the same
  if (nnz != other.nnz) {
    return false;
  }

  // check the matrix data against each other
  if (memcmp(vals, other.vals, nnz * sizeof(T)) != 0) {
    return false;
  }
  if (memcmp(innerIdx, other.innerIdx, nnz * sizeof(indexT)) != 0) {
    return false;
  }
  if (memcmp(outerPtr, other.outerPtr, (outerDim + 1) * sizeof(indexT)) != 0) {
    return false;
  }

  // if all the data is the same return true
  return true;
}

// Inequality Operator
template <typename T, typename indexT, bool columnMajor>
bool SparseMatrix<T, indexT, 1, columnMajor>::operator!=(const SparseMatrix<T, indexT, 1, columnMajor> &other) {
  return !(*this == other);
}

// Coefficent Access Operator
template <typename T, typename indexT, bool columnMajor>
T SparseMatrix<T, indexT, 1, columnMajor>::operator()(uint32_t row, uint32_t col) {
  
  #ifdef IVSPARSE_DEBUG
  // check if the row and column are in bounds
  if (row >= numRows || col >= numCols) {
    std::cerr << "Error: Index out of bounds" << std::endl;
    exit(1);
  }
  #endif

  // get the vector and index
  uint32_t vector = columnMajor ? col : row;
  uint32_t index = columnMajor ? row : col;

  // get an iterator for the desired vector
  for (typename SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator it(
           *this, vector);
       it; ++it) {
    if (it.getIndex() == (indexT)index) {
      // if the index is found return the value
      return it.value();
    }
  }

  // if the index is not found return 0
  return 0;
}

// Vector Access Operator
template <typename T, typename indexT, bool columnMajor>
typename SparseMatrix<T, indexT, 1, columnMajor>::Vector
SparseMatrix<T, indexT, 1, columnMajor>::operator[](uint32_t vec) {
  #ifdef IVSPARSE_DEBUG
  // check if the vector is out of bounds
  assert((vec < outerDim && vec >= 0) && "Vector index out of bounds");
  #endif

  // return a IVSparse vector
  typename IVSparse::SparseMatrix<T, indexT, 1, columnMajor>::Vector newVector(
      *this, vec);
  return newVector;
}

//* BLAS Operators *//

// Scalar Multiplication
template <typename T, typename indexT, bool columnMajor>
IVSparse::SparseMatrix<T, indexT, 1, columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::operator*(T scalar) {
  return scalarMultiply(scalar);
}

// In place scalar multiplication
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 1, columnMajor>::operator*=(T scalar) {
  return inPlaceScalarMultiply(scalar);
}

// IVSparse Matrix * IVSparse Vector Multiplication
template <typename T, typename indexT, bool columnMajor>
Eigen::VectorXd SparseMatrix<T, indexT, 1, columnMajor>::operator*(SparseMatrix<T, indexT, 1, columnMajor>::Vector &vec) {
  return vectorMultiply(vec);
}

// Matrix Vector Multiplication (IVSparse Eigen -> Eigen)
template <typename T, typename indexT, bool columnMajor>
Eigen::VectorXd SparseMatrix<T, indexT, 1, columnMajor>::operator*(Eigen::VectorXd &vec) {
  return vectorMultiply(vec);
}

// Matrix Matrix Multiplication (IVSparse Eigen -> Eigen)
template <typename T, typename indexT, bool columnMajor>
Eigen::Matrix<T, -1, -1> SparseMatrix<T, indexT, 1, columnMajor>::operator*(Eigen::Matrix<T, -1, -1> mat) {
  return matrixMultiply(mat);
}

}  // namespace IVSparse