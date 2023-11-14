/**
 * @file IVCSC_Vector_Methods.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief IVCSC Vector Methods
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

//* Constructors and Destructor *//

// Destructor
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::~Vector() {
  if (data != nullptr) {
    free(data);
  }
}

// Length Constructor
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::Vector(uint32_t length) {
  #ifdef IVSPARSE_DEBUG
    assert((length > 0) && "Vector length must be greater than 0");
  #endif

  this->length = length;
}

// IVSparse Matrix Constructor
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::Vector(
    IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor> &mat, uint32_t vec) {
  
  #ifdef IVSPARSE_DEBUG
  assert((vec >= 0 && vec < mat.outerSize()) && "Vector index out of bounds");
  assert((mat.outerSize() > 0 && mat.innerSize() > 0) && "Matrix is empty");
  #endif

  // get the length of the vector
  size = mat.getVectorSize(vec);
  length = mat.innerSize();

  // if the size is 0 then the vector is empty
  if (size == 0) {
    data = nullptr;
    endPtr = nullptr;
    return;
  }

  // set data pointer
  try {
    data = malloc(size);
  } catch (std::bad_alloc &e) {
    std::cerr << e.what() << '\n';
  }

  // copy the vector data into the vector
  memcpy(data, mat.vectorPointer(vec), size);

  // set the end pointer
  endPtr = (uint8_t *)data + size;

  // set the nnz
  if (nnz == 0 && size > 0) {
    // make an iterator for the vector
    IVSparse::SparseMatrix<T, indexT, compressionLevel,
                           columnMajor>::InnerIterator it(*this);

    // iterate through the vector until the index is found
    while (it) {
      nnz++;
      ++it;
    }
  }

}  // End of IVSparse Matrix Constructor

// Deep copy constructor
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::Vector(
    IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector &vec) {
  
  // set the size
  size = vec.size;

  // set the length
  length = vec.length;

  // if the size is 0 then the vector is empty
  if (size == 0) {
    data = nullptr;
    endPtr = nullptr;
    return;
  }

  // set data pointer
  try {
    data = malloc(size);
  } catch (std::bad_alloc &e) {
    std::cerr << e.what() << '\n';
  }

  // copy the vector data into the vector
  memcpy(data, vec.data, size);

  // set the end pointer
  endPtr = (uint8_t *)data + size;

  // set the nnz
  nnz = vec.nonZeros();

  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif
}

//* Private Class Methods *//

// User checks to confirm a valid vector
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
void SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::userChecks() {
  assert(std::is_floating_point<indexT>::value == false &&
         "The index type must be a non-floating point type");
  assert((compressionLevel == 1 || compressionLevel == 2 ||
          compressionLevel == 3) &&
         "The compression level must be either 1, 2, or 3");
  assert((std::is_arithmetic<T>::value && std::is_arithmetic<indexT>::value) &&
         "The value and index types must be numeric types");
  assert((std::is_same<indexT, bool>::value == false) &&
         "The index type must not be bool");
}

// Calculates the size of the vector in bytes
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
void SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::calculateCompSize() {
  size = 0;

  size += sizeof(void *);  // data pointer
  size += sizeof(void *);  // end pointer

  // add distance between data and end pointer in bytes to size
  size += (uint8_t *)endPtr - (uint8_t *)data;
}

//* Getters *//

// Get the inner dimension of the vector
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
uint32_t SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::innerSize() {
  return length;
}

// Get the outer dimension of the vector
template <typename T, typename indexT, uint8_t compressionLevel,bool columnMajor>
uint32_t SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::outerSize() {
  return 1;
}

// Get the number of non-zero elements in the vector
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
uint32_t SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::nonZeros() {
  return nnz;
}

// Get a pointer to the beginning of the vector
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
void *SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::begin() {
  return data;
}

// Get a pointer to the end of the vector
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
void *SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::end() {
  return endPtr;
}

// Get the size of the vector in bytes
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
size_t SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::byteSize() {
  return size;
}

// Update the value of the vector at the given index
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
T SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::coeff(uint32_t index) {
  
  #ifdef IVSPARSE_DEBUG
    assert(index < length && index >= 0 && "The index is out of bounds");
  #endif
  
  return (*this)[index];
}

// Get the length of the vector
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
uint32_t SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::getLength() {
  return length;
}

//* Utility Methods *//

// Prints the vector to console dense
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
void SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::print() {
  // if length is larger than 100 then print then don't print
  if (length > 100) {
    std::cout << "Vector is too large to print" << std::endl;
    return;
  }

  std::cout << "Vector: ";
  std::cout << std::endl;

  // print a dense vector
  for (uint32_t i = 0; i < length; i++) {
    std::cout << (*this)[i] << " ";
  }

  std::cout << std::endl;
}

//* Calculations *//

// Calculates the norm of the vector
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
inline double SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::norm() {
  double norm = 0;
  for (typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::InnerIterator it(*this);
       it; ++it) {
    norm += it.value() * it.value();
  }
  return sqrt(norm);
}

// Calculates the sum of the vector
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
inline T SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::sum() {
  double sum = 0;
  for (typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::InnerIterator it(*this);
       it; ++it) {
    sum += it.value();
  }
  return sum;
}

// Calculates the dot product of the vector with an Eigen dense vector
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
double SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::dot(Eigen::Matrix<T, -1, 1> &other) {
  
  double dot = 0;

  for (typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::InnerIterator it(*this);
       it; ++it) {
    dot += it.value() * other.coeff(it.row());
  }

  return dot;
}

// Calculates the dot product of the vector with an Eigen sparse vector
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
double SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::dot(Eigen::SparseVector<T, -1> &other) {
  
  double dot = 0;

  for (typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::InnerIterator it(*this);
       it; ++it) {
    dot += it.value() * other.coeff(it.row());
  }
  return dot;
}

//* Operator Overloads *//

// Assignment Operator
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector
SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::operator=(
    typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector &other) {
  
  // check if the vector is the same
  if (this == &other) {
    return *this;
  }

  // free the data if it is not null
  if (data != nullptr) {
    free(data);
  }

  // if the other vector is empty then return
  if (other.data == nullptr) {
    data = nullptr;
    endPtr = nullptr;
    size = 0;
    length = other.length;
    nnz = 0;
    return *this;
  }

  // copy the data
  data = (uint8_t *)malloc(other.size);
  memcpy(data, other.data, other.size);

  size = other.size;
  length = other.length;
  nnz = other.nnz;
  endPtr = (uint8_t *)data + size;

  // return this
  return *this;
}

// equality operator
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
bool SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::operator==(
    typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector &other) {
  
  // check if the length is the same
  if (length != other.length) {
    return false;
  }

  // check the nnz
  if (nnz != other.nnz) {
    return false;
  }

  // check data equality
  if (memcmp(data, other.data, size) != 0) {
    return false;
  }

  // return true
  return true;
}

// inequality operator
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
bool SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::operator!=( 
    typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector &other) {
  
  return !(*this == other);
}

// Coefficient Operator
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
T SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::operator[]( uint32_t index) {

  #ifdef IVSPARSE_DEBUG
  // check if the index is out of bounds
  assert(index < length && "The index is out of bounds");
  #endif

  // make an iterator for the vector
  IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::InnerIterator it(*this);

  // iterate through the vector until the index is found
  while (it) {
    if (it.getIndex() == (indexT)index) {
      return it.value();
    }
    ++it;
  }

  // if the index is not found then return 0
  return 0;
}

// In place scalar multiplication
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
void SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::operator*=(T scalar) {
  
  for (typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::InnerIterator it(*this); it; ++it) {
    if (it.isNewRun()) {
      it.coeff(it.value() * scalar);
    }
  }
}

// Scalar multiplication
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
typename IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector 
SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector::operator*(T scalar) {
  
  typename IVSparse::SparseMatrix<T, indexT, 2, columnMajor>::Vector newVector(*this);
  for (typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::InnerIterator it(newVector);
       it; ++it) {
    if (it.isNewRun()) {
      it.coeff(it.value() * scalar);
    }
  }
  return newVector;
}

}  // namespace IVSparse