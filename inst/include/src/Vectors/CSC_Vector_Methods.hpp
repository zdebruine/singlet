/**
 * @file CSC_Vector_Methods.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief CSC Vector Methods
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

//* Constructors and Destructor *//

// Destructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::Vector::~Vector() {
  if (vals != nullptr) {
    free(vals);
  }
  if (innerIdx != nullptr) {
    free(innerIdx);
  }
}

// IVSparse Matrix Constructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::Vector::Vector(
    IVSparse::SparseMatrix<T, indexT, 1, columnMajor> &mat, uint32_t vec) {

  #ifdef IVSPARSE_DEBUG
  // make sure the vector is in bounds
  assert((vec >= 0 && vec < mat.outerSize()) && "Vector index out of bounds");
  // make sure the matrix is not empty
  assert((mat.outerSize() > 0 && mat.innerSize() > 0) && "Matrix is empty");
  #endif

  length = mat.innerSize();

  // if the vector is empty, return
  if (mat.byteSize() == 0) {
    vals = nullptr;
    innerIdx = nullptr;
    return;
  }

  nnz = mat.outerPtr[vec + 1] - mat.outerPtr[vec];

  try {
    vals = (T *)malloc(nnz * sizeof(T));
    innerIdx = (indexT *)malloc(nnz * sizeof(indexT));
  } catch (std::bad_alloc &e) {
    std::cerr << "Allocation failed: " << e.what() << '\n';
  }

  memcpy(vals, mat.vals + mat.outerPtr[vec], nnz * sizeof(T));
  memcpy(innerIdx, mat.innerIdx + mat.outerPtr[vec], nnz * sizeof(indexT));

  calculateCompSize();

  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif
}

// IVSparse Vector Constructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::Vector::Vector(
    IVSparse::SparseMatrix<T, indexT, 1, columnMajor>::Vector &vec) {
  
  length = vec.length;
  nnz = vec.nnz;
  size = vec.size;

  try {
    vals = (T *)malloc(nnz * sizeof(T));
    innerIdx = (indexT *)malloc(nnz * sizeof(indexT));
  } catch (std::bad_alloc &e) {
    std::cerr << "Allocation failed: " << e.what() << '\n';
  }

  memcpy(vals, vec.vals, nnz * sizeof(T));
  memcpy(innerIdx, vec.innerIdx, nnz * sizeof(indexT));

  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif
}

//* Private Class Methods *//

// User Checks
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 1, columnMajor>::Vector::userChecks() {
  assert(std::is_floating_point<indexT>::value == false &&
         "The index type must be a non-floating point type");
  assert((std::is_arithmetic<T>::value && std::is_arithmetic<indexT>::value) &&
         "The value and index types must be numeric types");
  assert((std::is_same<indexT, bool>::value == false) &&
         "The index type must not be bool");
}

// Calculate the size of the vector
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 1, columnMajor>::Vector::calculateCompSize() {
  size = 0;

  // calculate the size of the vector
  size += sizeof(T) * nnz;
  size += sizeof(indexT) * nnz;
}

//* Getters *//

// Get the inner dimension of the vector
template <typename T, typename indexT, bool columnMajor>
uint32_t SparseMatrix<T, indexT, 1, columnMajor>::Vector::innerSize() {
  return length;
}

// Get the outer dimension of the vector
template <typename T, typename indexT, bool columnMajor>
uint32_t SparseMatrix<T, indexT, 1, columnMajor>::Vector::outerSize() {
  return 1;
}

// Get the number of nonzeros in the vector
template <typename T, typename indexT, bool columnMajor>
uint32_t SparseMatrix<T, indexT, 1, columnMajor>::Vector::nonZeros() {
  return nnz;
}

// Get the size of the vector in bytes
template <typename T, typename indexT, bool columnMajor>
size_t SparseMatrix<T, indexT, 1, columnMajor>::Vector::byteSize() {
  return size;
}

// Get a value from the vector at index
template <typename T, typename indexT, bool columnMajor>
T SparseMatrix<T, indexT, 1, columnMajor>::Vector::coeff(uint32_t index) {
  
  #ifdef IVSPARSE_DEBUG
    assert((index >= 0 && index < length) && "Index out of bounds");
  #endif
  
  return (*this)[index];
}

// Get the length of the vector
template <typename T, typename indexT, bool columnMajor>
uint32_t SparseMatrix<T, indexT, 1, columnMajor>::Vector::getLength() {
  return length;
}

// get values
template <typename T, typename indexT, bool columnMajor>
T *SparseMatrix<T, indexT, 1, columnMajor>::Vector::getValues() const {
  return vals;
}

// get inner indices
template <typename T, typename indexT, bool columnMajor>
indexT *SparseMatrix<T, indexT, 1, columnMajor>::Vector::getInnerIndices()
    const {
  return innerIdx;
}

//* Utility Methods *//

// print
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 1, columnMajor>::Vector::print() {
  std::cout << "Vector: " << std::endl;
  std::cout << std::endl;

  // print the denese vector up to 100 elements
  for (uint32_t i = 0; i < std::min(length, (uint32_t)100); i++) {
    std::cout << (*this)[i] << " ";
  }

  std::cout << std::endl;
}

//* Operator Overloads *//

// Assignment Operator
template <typename T, typename indexT, bool columnMajor>
typename SparseMatrix<T, indexT, 1, columnMajor>::Vector
SparseMatrix<T, indexT, 1, columnMajor>::Vector::operator=(
    typename SparseMatrix<T, indexT, 1, columnMajor>::Vector &vec) {
  
  // check if the vector is the same
  if (this == &vec) {
    return *this;
  }

  if (vals != nullptr) {
    delete[] vals;
  }
  if (innerIdx != nullptr) {
    delete[] innerIdx;
  }

  length = vec.length;
  nnz = vec.nnz;
  size = vec.size;

  // if the vector is empty, return
  if (size == 0) {
    vals = nullptr;
    innerIdx = nullptr;
    return *this;
  }

  try {
    vals = (T *)malloc(nnz * sizeof(T));
    innerIdx = (indexT *)malloc(nnz * sizeof(indexT));
  } catch (std::bad_alloc &e) {
    std::cerr << "Allocation failed: " << e.what() << '\n';
  }

  memcpy(vals, vec.vals, nnz * sizeof(T));
  memcpy(innerIdx, vec.innerIdx, nnz * sizeof(indexT));

  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif

  return *this;
}

// Equality Operator
template <typename T, typename indexT, bool columnMajor>
bool SparseMatrix<T, indexT, 1, columnMajor>::Vector::operator==(
    typename SparseMatrix<T, indexT, 1, columnMajor>::Vector &vec) {
  
  // check if the vector dimensions are the same
  if (length != vec.length || nnz != vec.nnz) {
    return false;
  }

  // check if the values and indices are the same
  for (uint32_t i = 0; i < nnz; i++) {
    if (vals[i] != vec.vals[i] || innerIdx[i] != vec.innerIdx[i]) {
      return false;
    }
  }

  // if all the values and indices are the same, return true
  return true;
}

// Inequality Operator
template <typename T, typename indexT, bool columnMajor>
bool SparseMatrix<T, indexT, 1, columnMajor>::Vector::operator!=(
    typename SparseMatrix<T, indexT, 1, columnMajor>::Vector &vec) {
  
  return !(*this == vec);
}

// Bracket Operator
template <typename T, typename indexT, bool columnMajor>
T SparseMatrix<T, indexT, 1, columnMajor>::Vector::operator[](uint32_t index) {
  #ifdef IVSPARSE_DEBUG
  assert((index >= 0 && index < length) && "Index out of bounds");
  #endif

  for (uint32_t i = 0; i < nnz; i++) {
    if (innerIdx[i] == index) {
      return vals[i];
    }
  }

  return 0;
}

}  // namespace IVSparse