/**
 * @file IVSparse_Base_Methods.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief IVSparse Sparse Matrix Base Methods
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    // Calculates the number of bytes needed to store a value
    inline uint8_t SparseMatrix::byteWidth(size_t size) {
        if (size <= 0xFF){
            return 1;
        }
        else if (size <= 0xFFFF){
            return 2;
        }
        else if (size <= 0xFFFFFF){
            return 3;
        }
        else if (size <= 0xFFFFFFFF){
            return 4;
        }
        else if (size <= 0xFFFFFFFFFF){
            return 5;
        }
        else if (size <= 0xFFFFFFFFFFFF){
            return 6;
        }
        else if (size <= 0xFFFFFFFFFFFFFF){
            return 7;
        }
        else{
            return 8;
        }

    }

    // Gets the number of rows in the matrix
    uint32_t SparseMatrix::rows() const { return numRows; }

    // Gets the number of columns in the matrix
    uint32_t SparseMatrix::cols() const { return numCols; }

    // Gets the inner dimension of the matrix
    uint32_t SparseMatrix::innerSize() const { return innerDim; }

    // Gets the outer dimension of the matrix
    uint32_t SparseMatrix::outerSize() const { return outerDim; }

    // Gets the number of non-zero elements in the matrix
    uint32_t SparseMatrix::nonZeros() const { return nnz; }

    // Gets the number of bytes needed to store the matrix
    size_t SparseMatrix::byteSize() const { return compSize; }

}  // namespace IVSparse