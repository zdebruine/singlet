/**
 * @file IVSparse_SparseMatrixBase.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief IVSparse Sparse Matrix Base Class Declerations
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    /**
     * IVSparse Sparse Matrix Base Class \n \n
     * This is the overarching base class for the IVSparse Sparse Matrix
     * Library. It contains methods and variables shared between all
     * compression levels of IVSparse Sparse Matrices and serves to reduce
     * code duplication.
     */
    class SparseMatrixBase {
        private:
        //* The Matrix Info *//

        uint32_t innerDim = 0;  // The inner dimension of the matrix
        uint32_t outerDim = 0;  // The outer dimension of the matrix

        uint32_t numRows = 0;  // The number of rows in the matrix
        uint32_t numCols = 0;  // The number of columns in the matrix

        uint32_t nnz = 0;  // The number of non-zero values in the matrix

        size_t compSize = 0;  // The size of the compressed matrix in bytes

        //* The Value and Index Types *//

        uint32_t val_t;  // Information about the value type (size, signededness, etc.)
        uint32_t index_t;  // Information about the index type (size)

        uint32_t* metadata = nullptr;  // The metadata of the matrix

        //* Private Methods *//

        // Calculates the number of bytes needed to store a value
        inline uint8_t byteWidth(size_t size);

        // Creates value type information
        virtual void encodeValueType() = 0;

        // Checks the value type information
        virtual void checkValueType() = 0;

        // User checks to confirm a valid matrix
        virtual void userChecks() = 0;

        // Calculates the size of the matrix in bytes
        virtual void calculateCompSize() = 0;

        public:
        //* Friends *//

        // IVSparse Sparse Matrix Class
        template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
        friend class SparseMatrix;

        //* Constructors *//

        // Default Constructor
        SparseMatrixBase() {};

        //* Getters *//

        /**
         * @returns The number of rows in the matrix.
         */
        uint32_t rows() const;

        /**
         * @returns The number of columns in the matrix.
         */
        uint32_t cols() const;

        /**
         * @returns The inner dimension of the matrix.
         */
        uint32_t innerSize() const;

        /**
         * @returns The outer dimension of the matrix.
         */
        uint32_t outerSize() const;

        /**
         * @returns The number of non-zero elements in the matrix.
         */
        uint32_t nonZeros() const;

        /**
         * @returns The size of the matrix in bytes.
         */
        uint64_t byteSize() const;

        //* Utility Methods *//

        /**
         * Writes the matrix to a file with the given filename.
         */
        virtual void write(const char* filename) = 0;

        /**
         * Prints the matrix to the console.
         */
        virtual void print() = 0;

    };  // class SparseMatrixBase

}  // namespace IVSparse