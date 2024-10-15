/**
 * @file VCSC_Vector_Methods.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Vector Methods for VCSC Sparse Matrices
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

//* Constructors and Destructor *//

// Destructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 2, columnMajor>::Vector::~Vector() {
}

// IVSparse Matrix Constructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 2, columnMajor>::Vector::Vector(
    IVSparse::SparseMatrix<T, indexT, 2, columnMajor> &mat, uint32_t vec) {
  
  #ifdef IVSPARSE_DEBUG
  // make sure the vector is in bounds
  assert((vec >= 0 && vec < mat.outerSize()) && "Vector index out of bounds");

  // make sure the matrix is not empty
  assert((mat.outerSize() > 0 && mat.innerSize() > 0) && "Matrix is empty");
  #endif

  // check if the vector is empty
  if (mat.getNumUniqueVals(vec) == 0) {
    size = 0;
    length = mat.innerSize();
    return;
  }

  length = mat.innerSize();
  nnz = mat.getNumIndices(vec);


  // copy the data[vec] to the vector
  data = mat.getMap(vec);

  calculateCompSize();
}

// Deep copy constructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 2, columnMajor>::Vector::Vector(
    IVSparse::SparseMatrix<T, indexT, 2, columnMajor>::Vector &vec) {
  
  // set the variables
  length = vec.length;
  size = vec.size;
  nnz = vec.nnz;
  indexWidth = vec.indexWidth;

  if (size == 0) {
    return;
  }

  // copy the data
  data = vec.data;

  // user checks
  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif
}

//* Private Class Methods *//

// User checks to confirm a valid vector
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 2, columnMajor>::Vector::userChecks() {
  assert(std::is_floating_point<indexT>::value == false &&
         "The index type must be a non-floating point type");
  assert((std::is_arithmetic<T>::value && std::is_arithmetic<indexT>::value) &&
         "The value and index types must be numeric types");
  assert((std::is_same<indexT, bool>::value == false) &&
         "The index type must not be bool");
}

// Calculates the size of the vector in bytes
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 2, columnMajor>::Vector::calculateCompSize() {
  size = 0;

  size += sizeof(T) * data.size();       // values
  size += sizeof(indexT) * data.size();  // counts
  for (auto j = 0; j < data.size(); j++) {
    size += sizeof(indexT) * data[j].size();  // indices
  }
  size += sizeof(indexT);  // len
}

//* Getters *//

// Get the inner size of the vector
template <typename T, typename indexT, bool columnMajor>
uint32_t SparseMatrix<T, indexT, 2, columnMajor>::Vector::innerSize() {
  return length;
}

// Get the outer size of the vector
template <typename T, typename indexT, bool columnMajor>
uint32_t SparseMatrix<T, indexT, 2, columnMajor>::Vector::outerSize() {
  return 1;
}

// Get the number of non-zero elements in the vector
template <typename T, typename indexT, bool columnMajor>
uint32_t SparseMatrix<T, indexT, 2, columnMajor>::Vector::nonZeros() {
  return nnz;
}

// Get the byte size of the vector
template <typename T, typename indexT, bool columnMajor>
size_t SparseMatrix<T, indexT, 2, columnMajor>::Vector::byteSize() {
  return size;
}

// Get the value at the given index
template <typename T, typename indexT, bool columnMajor>
T SparseMatrix<T, indexT, 2, columnMajor>::Vector::coeff(uint32_t index) {
  
  #ifdef IVSPARSE_DEBUG
    // check if the index is out of bounds
    assert(index < length && index >= 0 && "The index is out of bounds");
  #endif
  
  return (*this)[index];
}

// Get the length of the vector
template <typename T, typename indexT, bool columnMajor>
uint32_t SparseMatrix<T, indexT, 2, columnMajor>::Vector::getLength() {
  return length;
}

// Get a pointer to the values of the vector
template <typename T, typename indexT, bool columnMajor>
std::vector<T> SparseMatrix<T, indexT, 2, columnMajor>::Vector::getValues() {
  std::vector<T> temp(length);

  for (auto& i : data) {
    temp.push_back(i.value);
  }

  return temp;
}

// Get a pointer to the counts of the vector
template <typename T, typename indexT, bool columnMajor>
std::vector<indexT> SparseMatrix<T, indexT, 2, columnMajor>::Vector::getCounts() {
  std::vector<indexT> temp(length);

  for (auto& [key, value] : data) {
    temp.push_back(value.size());
  }

  return temp;
}

// Get a pointer to the indices of the vector
template <typename T, typename indexT, bool columnMajor>
std::vector<indexT> SparseMatrix<T, indexT, 2, columnMajor>::Vector::getIndices() {
  std::vector<indexT> temp(nnz);

  for (auto& [key, value] : data) {
    for (auto& i : value) {
      temp.push_back(i);
    }
  }

  return temp;
}

// Get the number of unique values in the vector
template <typename T, typename indexT, bool columnMajor>
indexT SparseMatrix<T, indexT, 2, columnMajor>::Vector::uniqueVals() {
  return data.size();
}

// returns the innner map
template <typename T, typename indexT, bool columnMajor>
std::map<T, std::vector<indexT>> SparseMatrix<T, indexT, 2, columnMajor>::Vector::getData() {
  return data;
}

//* Utility Methods *//

// Print the vector to console
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 2, columnMajor>::Vector::print() {
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

//* Calculation Methods *//

// Calculate the norm of the vector
template <typename T, typename indexT, bool columnMajor>
inline double SparseMatrix<T, indexT, 2, columnMajor>::Vector::norm() {
  double norm = 0;
  for (typename SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator it(*this); it; ++it) {
    norm += it.value() * it.value();
  }
  return sqrt(norm);
}

// Calculate the sum of the vector
template <typename T, typename indexT, bool columnMajor>
inline T SparseMatrix<T, indexT, 2, columnMajor>::Vector::sum() {
  double sum = 0;
  for (typename SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator it(*this); it; ++it) {
    sum += it.value();
  }
  return sum;
}

// Calculate the dot product with an Eigen dense vector
template <typename T, typename indexT, bool columnMajor>
double SparseMatrix<T, indexT, 2, columnMajor>::Vector::dot(Eigen::Matrix<T, -1, 1> &other) {
  
  double dot = 0;

  for (typename SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator it(
           *this);
       it; ++it) {
    dot += it.value() * other.coeff(it.row());
  }

  return dot;
}

// Calculate the dot product with an Eigen sparse vector
template <typename T, typename indexT, bool columnMajor>
double SparseMatrix<T, indexT, 2, columnMajor>::Vector::dot(Eigen::SparseVector<T, -1> &other) {
  double dot = 0;

  for (typename SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator it(
           *this);
       it; ++it) {
    dot += it.value() * other.coeff(it.row());
  }

  return dot;
}

//* Operator Overloads *//

template <typename T, typename indexT, bool columnMajor>
typename SparseMatrix<T, indexT, 2, columnMajor>::Vector
SparseMatrix<T, indexT, 2, columnMajor>::Vector::operator=(typename SparseMatrix<T, indexT, 2, columnMajor>::Vector &other) {
  
  // check if the vector is the same
  if (this == &other) {
    return *this;
  }

  // set the variables
  length = other.length;
  size = other.size;
  nnz = other.nnz;
  indexWidth = other.indexWidth;

  if (size == 0) {
    return *this;
  }

  // copy the data
  data = other.data;

  // user checks
  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif

  // return this
  return *this;
}

// equality operator
template <typename T, typename indexT, bool columnMajor>
bool SparseMatrix<T, indexT, 2, columnMajor>::Vector::operator==(
    typename SparseMatrix<T, indexT, 2, columnMajor>::Vector &other) {
  
  // check if the lengths are the same
  if (length != other.length) {
    return false;
  }

  // check if the nnz are the same
  if (nnz != other.nnz) {
    return false;
  }

  // check if the values are the same
  if (data != other.data) {
    return false;
  }

  // return true if all checks pass
  return true;
}

// inequality operator
template <typename T, typename indexT, bool columnMajor>
bool SparseMatrix<T, indexT, 2, columnMajor>::Vector::operator!=(
    typename SparseMatrix<T, indexT, 2, columnMajor>::Vector &other) {
  return !(*this == other);
}

// coefficient access operator
template <typename T, typename indexT, bool columnMajor>
T SparseMatrix<T, indexT, 2, columnMajor>::Vector::operator[](uint32_t index) {
  
  #ifdef IVSPARSE_DEBUG
  // check if the index is out of bounds
  assert(index < length && "The index is out of bounds");
  #endif

  // make an iterator for the vector
  IVSparse::SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator it(*this);

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

// Scalar multiplication operator
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 2, columnMajor>::Vector::operator*=(T scalar) {
  std::map<T, std::vector<indexT>> newValues;
  for (const auto& [key, value] : data) {
    newValues[key * scalar] = value;
  }
  data = newValues;
}

// Scalar multiplication operator (returns a new vector)
template <typename T, typename indexT, bool columnMajor>
typename IVSparse::SparseMatrix<T, indexT, 2, columnMajor>::Vector
SparseMatrix<T, indexT, 2, columnMajor>::Vector::operator*(T scalar) {
  
  typename IVSparse::SparseMatrix<T, indexT, 2, columnMajor>::Vector newVector(*this);

  std::map<T, std::vector<indexT>> newValues;
  for (const auto &[key, value] : newVector.data) {
    newValues[key * scalar] = value;
  }
  newVector.data = newValues;

  return newVector;
}

}  // namespace IVSparse