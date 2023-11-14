/**
 * @file CSC_Private_Methods.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Private Methods for CSC Sparse Matrices
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

// Calculates the number of bytes needed to store a value
template <typename T, typename indexType, bool columnMajor>
inline uint8_t SparseMatrix<T, indexType, 1, columnMajor>::byteWidth(size_t size) {
    if (size <= 0xFF) {
        return 1;
    } else if (size <= 0xFFFF) {
        return 2;
    } else if (size <= 0xFFFFFF) {
        return 3;
    } else if (size <= 0xFFFFFFFF) {
        return 4;
    } else if (size <= 0xFFFFFFFFFF) {
        return 5;
    } else if (size <= 0xFFFFFFFFFFFF) {
        return 6;
    } else if (size <= 0xFFFFFFFFFFFFFF) {
        return 7;
    } else {
        return 8;
    }
}

// Encodes the value type of the matrix in a uint32_t
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 1, columnMajor>::encodeValueType() {
    uint8_t byte0 = sizeof(T);
    uint8_t byte1 = std::is_floating_point<T>::value ? 1 : 0;
    uint8_t byte2 = std::is_signed<T>::value ? 1 : 0;
    uint8_t byte3 = columnMajor ? 1 : 0;

    val_t = (byte3 << 24) | (byte2 << 16) | (byte1 << 8) | byte0;
}

// Checks if the value type is correct for the matrix
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 1, columnMajor>::checkValueType() {
    uint8_t byte0 = val_t & 0xFF;
    uint8_t byte1 = (val_t >> 8) & 0xFF;
    uint8_t byte2 = (val_t >> 16) & 0xFF;
    uint8_t byte3 = (val_t >> 24) & 0xFF;
    assert(byte0 == sizeof(T) && "Value type size does not match");
    assert(byte1 == std::is_floating_point<T>::value &&
           "Value type is not floating point");
    assert(byte2 == std::is_signed<T>::value && "Value type is not signed");
    assert(byte3 == columnMajor && "Major direction does not match");
}

// performs some simple user checks on the matrices metadata
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 1, columnMajor>::userChecks() {
    assert((innerDim > 1 || outerDim > 1 || nnz > 1) &&
           "The matrix must have at least one row, column, and nonzero value");
    assert(std::is_floating_point<indexT>::value == false &&
           "The index type must be a non-floating point type");
    assert((std::is_arithmetic<T>::value && std::is_arithmetic<indexT>::value) &&
           "The value and index types must be numeric types");
    assert((std::is_same<indexT, bool>::value == false) &&
           "The index type must not be bool");
    assert((innerDim < std::numeric_limits<indexT>::max() &&
            outerDim < std::numeric_limits<indexT>::max()) &&
           "The number of rows and columns must be less than the maximum value "
           "of the index type");
    checkValueType();
}

// Calculates the current byte size of the matrix in memory
template <typename T, typename indexT, bool columnMajor>
void SparseMatrix<T, indexT, 1, columnMajor>::calculateCompSize() {
    // set compSize to zero
    compSize = 0;

    // add the size of the metadata
    compSize += META_DATA_SIZE;

    // add the csc vectors
    compSize += sizeof(T) * nnz;                  // values
    compSize += sizeof(indexT) * nnz;             // innerIdx
    compSize += sizeof(indexT) * (outerDim + 1);  // outerPtr
}

}  // namespace IVSparse