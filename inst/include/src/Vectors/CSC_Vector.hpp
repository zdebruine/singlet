/**
 * @file CSC_Vector.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief CSC Vector Class Declerations
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

/**
 * CSC Vector Class \n \n
 * The CSC Vector class is a vector class that is used to work with
 * CSC matrices. It works with the same logic as the corresponding
 * matrix compression level and is useful when working with these matrices.
 */
template <typename T, typename indexT, bool columnMajor>
class SparseMatrix<T, indexT, 1, columnMajor>::Vector {
 private:
  //* Private Class Variables *//

  size_t size = 0;  // size of the vector in bytes

  T *vals = nullptr;           // values of the vector
  indexT *innerIdx = nullptr;  // inner indices of the vector

  uint32_t length = 0;  // length of the vector
  uint32_t nnz = 0;     // number of non-zero elements in the vector

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
   * IVSparse Matrix to Vector Constructor \n \n
   * Creates a vector from a CSC Matrix at the given vector index.
   *
   * @note Can only get a vector from a matrix in the storage order of the
   * matrix.
   */
  Vector(IVSparse::SparseMatrix<T, indexT, 1, columnMajor> &mat, uint32_t vec);

  /**
   * Deep Copy Vector Constructor \n \n
   * Creates a deep copy of the given vector.
   */
  Vector(IVSparse::SparseMatrix<T, indexT, 1, columnMajor>::Vector &vec);

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

  /**
   * @returns A pointer to the values of the vector.
   */
  T *getValues() const;

  /**
   * @returns A pointer to the inner indices of the vector.
   */
  indexT *getInnerIndices() const;

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

  //* Operator Overloads *//

  // Coefficient Access Operator
  T operator[](uint32_t index);

  // Assignment Operator
  typename SparseMatrix<T, indexT, 1, columnMajor>::Vector operator=(
      typename SparseMatrix<T, indexT, 1, columnMajor>::Vector &vec);

  // Equality Operator
  bool operator==(
      typename SparseMatrix<T, indexT, 1, columnMajor>::Vector &vec);

  // Inequality Operator
  bool operator!=(
      typename SparseMatrix<T, indexT, 1, columnMajor>::Vector &vec);

};  // class Vector

}  // namespace IVSparse