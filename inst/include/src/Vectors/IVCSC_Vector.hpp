/**
 * @file IVCSC_Vector.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief IVCSC Vector Class Declerations
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

/**
 * @tparam T Type of the values in the matrix
 * @tparam indexT Type of the indices in the matrix
 * @tparam compressionLevel Compression level of the matrix
 * @tparam columnMajor Storage order of the matrix
 *
 * IVCSC Vector Class \n \n
 * The IVCSC Vector class is a vector class that is used to work with
 * IVCSC matrices. It works with the same logic as the corresponding
 * matrix compression level and is useful when working with these matrices.
 */
template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
class SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector {
 private:
  //* Private Class Variables *//

  size_t size = 0;  // size of the vector in bytes

  void *data = nullptr;    // data of the vector
  void *endPtr = nullptr;  // pointer to the end of the vector

  uint32_t length = 0;  // length of the vector

  uint8_t indexWidth = 1;  // width of the indices

  uint32_t nnz = 0;  // number of non-zero elements in the vector

  //* Private Class Methods *//

  // User checks to confirm a valid vector
  void userChecks();

  // Calculates the size of the vector in bytes
  void calculateCompSize();

 public:
  //* Constructors & Destructor *//
  /** @name Constructors
   */
  ///@{

  /**
   * Default Vector Constructor \n \n
   * Creates an empty vector with everything set to null/zero.
   */
  Vector(){};

  /**
   * Length Vector Constructor \n \n
   * Creates a vector of the given length with everything set to null/zero.
   */
  Vector(uint32_t length);

  /**
   * IVSparse Matrix to Vector Constructor \n \n
   * Creates a vector from a IVCSC Matrix at the given vector index.
   *
   * @note Can only get a vector from a matrix in the storage order of the
   * matrix.
   */
  Vector(IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor> &mat, uint32_t vec);

  /**
   * Deep Copy Vector Constructor \n \n
   * Creates a deep copy of the given vector.
   */
  Vector(IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector &vec);

  /**
   * Destroys the vector.
   */
  ~Vector();

  ///@}

  //* Getters *//
  /** @name Getters
   */
  ///@{

  /**
   * @returns The coefficient at the given index.
   */
  T coeff(uint32_t index);

  /**
   * @returns A pointer to the beginning of the vector.
   */
  void *begin();

  /**
   * @returns A pointer to the end of the vector.
   */
  void *end();

  /**
   * @returns The size of the vector in bytes.
   */
  size_t byteSize();

  /**
   * @returns The inner size of the vector.
   */
  uint32_t innerSize();

  /**
   * @returns The outer size of the vector.
   */
  uint32_t outerSize();

  /**
   * @returns The number of non-zero elements in the vector.
   */
  uint32_t nonZeros();

  /**
   * @returns The length of the vector.
   */
  uint32_t getLength();

  ///@}

  //* Utility Methods *//
  /** @name Utility Methods
   */
  ///@{

  /**
   * Prints the vector dense to the console.
   */
  void print();

  ///@}

  //* Calculations *//
  /** @name Calculation Methods
   */
  ///@{

  /**
   * @returns The norm of the vector.
   */
  double norm();

  /**
   * @returns The sum of the vector.
   */
  T sum();

  /**
   * @returns The dot product of the vector and an Eigen Dense Vector.
   */
  double dot(Eigen::Matrix<T, -1, 1> &other);

  /**
   * @returns The dot product of the vector and an Eigen Sparse Vector.
   */
  double dot(Eigen::SparseVector<T, -1> &other);

  ///@}

  //* Operator Overloads *//

  // In place scalar multiplication
  void operator*=(T scalar);

  // scalar multiplication
  typename IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector operator*(T scalar);

  // equality operator
  bool operator==(typename SparseMatrix<T, indexT, compressionLevel,
                                        columnMajor>::Vector &vec);

  // inequality operator
  bool operator!=(typename SparseMatrix<T, indexT, compressionLevel,
                                        columnMajor>::Vector &vec);

  // coefficient access
  T operator[](uint32_t index);

  // boolean operator
  operator bool() { return (char *)endPtr - indexWidth > data; };

  // assignment operator
  typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector
  operator=(typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector &vec);

};  // class Vector

}  // namespace IVSparse