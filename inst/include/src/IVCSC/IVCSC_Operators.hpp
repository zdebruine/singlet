/**
 * @file IVCSC_Operators.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Operator Overloads for IVCSC Sparse Matrices
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    // Assignment Operator
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>& SparseMatrix<T, indexT, compressionLevel, columnMajor>::operator=(const IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>& other) {

        if (this != &other) {
            // free old data
            if (data != nullptr) {
                for (uint32_t i = 0; i < outerDim; i++) {
                    if (data[i] != nullptr) {
                        free(data[i]);
                    }
                }
                free(data);
            }
            if (endPointers != nullptr) {
                free(endPointers);
            }
            if (metadata != nullptr) {
                delete[] metadata;
            }

            // set the dimensions
            numRows = other.numRows;
            numCols = other.numCols;
            outerDim = other.outerDim;
            innerDim = other.innerDim;
            nnz = other.nnz;
            compSize = other.compSize;

            // allocate the memory
            try {
                data = (void**)malloc(outerDim * sizeof(void*));
                endPointers = (void**)malloc(outerDim * sizeof(void*));
                metadata = new uint32_t[NUM_META_DATA];
            }
            catch (std::bad_alloc& e) {
                std::cerr << "Error: Could not allocate memory for IVSparse matrix"
                    << std::endl;
                exit(1);
            }

            // copy the metadata
            memcpy(metadata, other.metadata, sizeof(uint32_t) * NUM_META_DATA);

            // set the index and value types
            encodeValueType();
            index_t = other.index_t;

            // copy the data
            for (uint32_t i = 0; i < outerDim; i++) {
                // if the vector is empty, set the data pointer to nullptr
                if (other.data[i] == nullptr) {
                    data[i] = nullptr;
                    endPointers[i] = nullptr;
                    continue;
                }

                try {
                    data[i] = malloc(other.getVectorSize(i));
                }
                catch (std::bad_alloc& e) {
                    std::cerr << "Error: Could not allocate memory for IVSparse matrix"
                        << std::endl;
                    exit(1);
                }

                memcpy(data[i], other.data[i], other.getVectorSize(i));
                endPointers[i] = (uint8_t*)data[i] + other.getVectorSize(i);
            }
        }
        return *this;
    }

    // Equality Operator
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    bool SparseMatrix<T, indexT, compressionLevel, columnMajor>::operator==(const SparseMatrix<T, indexT, compressionLevel, columnMajor>& other) const {
        // bool SparseMatrix<T, indexT, 2               , columnMajor>::operator==(const SparseMatrix<T, indexT,                2, columnMajor>& other) const {

            // check if the two matrices are equal

            // first check the metadata using memcompare
        if (memcmp(metadata, other.metadata, sizeof(uint32_t) * NUM_META_DATA) != 0)
            return false;

        // iterate through the data and compare each element
        for (uint32_t i = 0; i < outerDim; i++) {
            if (memcmp(data[i], other.data[i], getVectorSize(i)) != 0) return false;
        }

        return true;
    }

    // Inequality Operator
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    bool SparseMatrix<T, indexT, compressionLevel, columnMajor>::operator!=(const SparseMatrix<T, indexT, compressionLevel, columnMajor>& other) {

        return !(*this == other);
    }

    // Coefficent Access Operator
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    T SparseMatrix<T, indexT, compressionLevel, columnMajor>::operator()(uint32_t row, uint32_t col) {
        #ifdef IVSPARSE_DEBUG
        // check if the row and column are in bounds
        if (row >= numRows || col >= numCols) {
            std::cerr << "Error: Index out of bounds" << std::endl;
            exit(1);
        }
        #endif

        uint32_t vector = columnMajor ? col : row;
        uint32_t index = columnMajor ? row : col;

        // if the vector is empty return 0
        if (data[vector] == nullptr) return 0;

        // get an iterator for the desired vector
        for (typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::InnerIterator it(*this, vector); it; ++it) {
            if (it.getIndex() == (indexT)index) {
                return it.value();
            }
        }

        // if the index is not found return 0
        return 0;
    }

    // Vector Access Operator
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector
        SparseMatrix<T, indexT, compressionLevel, columnMajor>::operator[](uint32_t vec) {

        #ifdef IVSPARSE_DEBUG
        // check if the vector is out of bounds
        assert((vec < outerDim && vec >= 0) && "Vector index out of bounds");
        #endif

        // return a IVSparse vector
        IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector newVector(*this, vec);

        return newVector;
    }

    //* BLAS Operators *//

    // Scalar Multiplication
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor> SparseMatrix<T, indexT, compressionLevel, columnMajor>::operator*(T scalar) {
        return scalarMultiply(scalar);
    }

    // In place scalar multiplication
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    void SparseMatrix<T, indexT, compressionLevel, columnMajor>::operator*=(T scalar) {
        return inPlaceScalarMultiply(scalar);
    }

    // IVSparse Matrix * IVSparse Vector Multiplication
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    Eigen::Matrix<T, -1, 1> SparseMatrix<T, indexT, compressionLevel, columnMajor>::operator*(
        SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector& vec) {

        return vectorMultiply(vec);
    }

    // Matrix Vector Multiplication (IVSparse Eigen -> Eigen)
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    Eigen::Matrix<T, -1, 1> SparseMatrix<T, indexT, compressionLevel, columnMajor>::operator*(
        Eigen::Matrix<T, -1, 1>& vec) {

        return vectorMultiply(vec);
    }

    // Matrix Matrix Multiplication (IVSparse Eigen -> Eigen)
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    Eigen::Matrix<T, -1, -1> SparseMatrix<T, indexT, compressionLevel, columnMajor>::operator*(Eigen::Matrix<T, -1, -1>& mat) {

        #ifdef IVSPARSE_DEBUG
        // check that the matrix is the correct size
        if (mat.rows() != numCols)
            throw std::invalid_argument(
                "The left matrix must have the same # of rows as columns in the right "
                "matrix!");
        #endif

        Eigen::Matrix<T, -1, -1> newMatrix = Eigen::Matrix<T, -1, -1>::Zero(mat.cols(), numRows);
        Eigen::Matrix<T, -1, -1> matTranspose = mat.transpose();

        // #ifdef IVSPARSE_HAS_OPENMP
        // #pragma omp parallel for
        // #endif
        for (uint32_t col = 0; col < numCols; col++) {
            for (typename SparseMatrix<T, indexT, 3, columnMajor>::InnerIterator matIter(*this, col); matIter; ++matIter) {
                newMatrix.col(matIter.row()) += matTranspose.col(col) * matIter.value();
            }
        }
        return newMatrix.transpose();
    }

} // namespace IVSparse