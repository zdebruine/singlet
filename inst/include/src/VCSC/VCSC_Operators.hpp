/**
 * @file VCSC_Operators.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Operator Overloads for VCSC Sparse Matrices
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    // Assignment Operator
    template <typename T, typename indexT, bool columnMajor>
    SparseMatrix<T, indexT, 2, columnMajor>& SparseMatrix<T, indexT, 2, columnMajor>::operator=(const IVSparse::SparseMatrix<T, indexT, 2, columnMajor>& other) {
        // check if the matrices are the same
        if (this != &other) {
            // free the old data
            if (values != nullptr) {
                for (uint32_t i = 0; i < outerDim; i++) {
                    if (values[i] != nullptr) free(values[i]);
                }
                free(values);
            }
            if (counts != nullptr) {
                for (uint32_t i = 0; i < outerDim; i++) {
                    if (counts[i] != nullptr) free(counts[i]);
                }
                free(counts);
            }
            if (indices != nullptr) {
                for (uint32_t i = 0; i < outerDim; i++) {
                    if (indices[i] != nullptr) free(indices[i]);
                }
                free(indices);
            }
            // free the metadata and size arrays
            if (valueSizes != nullptr) {
                free(valueSizes);
            }
            if (indexSizes != nullptr) {
                free(indexSizes);
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
                values = (T**)malloc(outerDim * sizeof(T*));
                counts = (indexT**)malloc(outerDim * sizeof(indexT*));
                indices = (indexT**)malloc(outerDim * sizeof(indexT*));

                valueSizes = (indexT*)malloc(outerDim * sizeof(indexT));
                indexSizes = (indexT*)malloc(outerDim * sizeof(indexT));

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
                try {
                    values[i] = (T*)malloc(other.valueSizes[i] * sizeof(T));
                    counts[i] = (indexT*)malloc(other.valueSizes[i] * sizeof(indexT));
                    indices[i] = (indexT*)malloc(other.indexSizes[i] * sizeof(indexT));
                }
                catch (std::bad_alloc& e) {
                    std::cerr << "Error: Could not allocate memory for IVSparse matrix"
                        << std::endl;
                    exit(1);
                }

                // copy the data
                memcpy(values[i], other.values[i], sizeof(T) * other.valueSizes[i]);
                memcpy(counts[i], other.counts[i], sizeof(indexT) * other.valueSizes[i]);
                memcpy(indices[i], other.indices[i],
                       sizeof(indexT) * other.indexSizes[i]);

                valueSizes[i] = other.valueSizes[i];
                indexSizes[i] = other.indexSizes[i];
            }
        }

        // return the new matrix
        return *this;

    }  // end assignment operator

    // Equality Operator
    template <typename T, typename indexT, bool columnMajor>
    bool SparseMatrix<T, indexT, 2, columnMajor>::operator==(const SparseMatrix<T, indexT, 2, columnMajor>& other) const {
        // first check the metadata using memcompare
        if (memcmp(metadata, other.metadata, sizeof(uint32_t) * NUM_META_DATA) != 0) {
            return false;
        }

        // check the value array
        for (uint32_t i = 0; i < outerDim; i++) {
            if (memcmp(values[i], other.values[i], sizeof(T) * valueSizes[i]) != 0) {
                return false;
            }
        }

        // check the index array
        for (uint32_t i = 0; i < outerDim; i++) {
            if (memcmp(indices[i], other.indices[i], sizeof(indexT) * indexSizes[i]) !=
                0) {
                return false;
            }
        }

        // check the count array
        for (uint32_t i = 0; i < outerDim; i++) {
            if (memcmp(counts[i], other.counts[i], sizeof(indexT) * valueSizes[i]) !=
                0) {
                return false;
            }
        }

        // if all of the above checks pass then the matrices are equal
        return true;
    }

    // Inequality Operator
    template <typename T, typename indexT, bool columnMajor>
    bool SparseMatrix<T, indexT, 2, columnMajor>::operator!=(const SparseMatrix<T, indexT, 2, columnMajor>& other) {
        return !(*this == other);
    }

    // Coefficent Access Operator
    template <typename T, typename indexT, bool columnMajor>
    T SparseMatrix<T, indexT, 2, columnMajor>::operator()(uint32_t row, uint32_t col) {
        #ifdef IVSPARSE_DEBUG
        // check if the row and column are in bounds
        assert((row < numRows && row >= 0) && "Row index out of bounds");
        assert((col < numCols && col >= 0) && "Column index out of bounds");
        #endif

        #ifdef IVSPARSE_DEBUG
        // check if the row and column are in bounds
        if (row >= numRows || col >= numCols) {
            std::cerr << "Error: Index out of bounds" << std::endl;
            exit(1);
        }
        #endif

        uint32_t vector = columnMajor ? col : row;
        uint32_t index = columnMajor ? row : col;

        // get an iterator for the desired vector
        for (typename SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator it(
            *this, vector);
            it; ++it) {
            if (it.getIndex() == (indexT)index) {
                return it.value();
            }
        }

        // if the index is not found return 0
        return 0;
    }

    // Vector Access Operator
    template <typename T, typename indexT, bool columnMajor>
    typename SparseMatrix<T, indexT, 2, columnMajor>::Vector SparseMatrix<T, indexT, 2, columnMajor>::operator[](uint32_t vec) {
        #ifdef IVSPARSE_DEBUG
        // check if the vector is out of bounds
        assert((vec < outerDim && vec >= 0) && "Vector index out of bounds");
        #endif

        // return a IVSparse vector
        typename IVSparse::SparseMatrix<T, indexT, 2, columnMajor>::Vector newVector(
            *this, vec);

        return newVector;
    }


    //* BLAS Operators *//

    // Scalar Multiplication
    template <typename T, typename indexT, bool columnMajor>
    IVSparse::SparseMatrix<T, indexT, 2, columnMajor> SparseMatrix<T, indexT, 2, columnMajor>::operator*(T scalar) {
        return scalarMultiply(scalar);
    }

    // In place scalar multiplication
    template <typename T, typename indexT, bool columnMajor>
    void SparseMatrix<T, indexT, 2, columnMajor>::operator*=(T scalar) {
        return inPlaceScalarMultiply(scalar);
    }

    // // IVSparse Matrix * IVSparse Vector Multiplication
    // template <typename T, typename indexT, bool columnMajor>
    // Eigen::Matrix<T, -1, 1>SparseMatrix<T, indexT, 2,
    // columnMajor>::operator*(SparseMatrix<T, indexT, 2, columnMajor>::Vector& vec)
    // {
    //     return vectorMultiply(vec);
    // }

    // Matrix Vector Multiplication (IVSparse Eigen -> Eigen)
    template <typename T, typename indexT, bool columnMajor>
    Eigen::Matrix<T, -1, 1> SparseMatrix<T, indexT, 2, columnMajor>::operator*(Eigen::Matrix<T, -1, 1>& vec) {
        return vectorMultiply(vec);
    }

    // Matrix Matrix Multiplication (IVSparse Eigen -> Eigen)
    template <typename T, typename indexT, bool columnMajor>
    Eigen::Matrix<T, -1, -1> SparseMatrix<T, indexT, 2, columnMajor>::operator* (Eigen::Matrix<T, -1, -1>& mat) {
        #ifdef IVSPARSE_DEBUG
        // check that the matrix is the correct size
        if (mat.rows() != numCols)
            throw std::invalid_argument(
                "The left matrix must have the same # of rows as columns in the right "
                "matrix!");
        #endif

        Eigen::Matrix<T, -1, -1> newMatrix = Eigen::Matrix<T, -1, -1>::Zero(mat.cols(), numRows);
        Eigen::Matrix<T, -1, -1> matTranspose = mat.transpose();

        // Fix Parallelism issue (race condition because of partial sums and
        // orientation of Sparse * Dense)
        for (uint32_t col = 0; col < numCols; col++) {
            for (typename SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator
                 matIter(*this, col);
                 matIter; ++matIter) {
                newMatrix.col(matIter.row()) += matTranspose.col(col) * matIter.value();
            }
        }
        return newMatrix.transpose();
    }

}  // end namespace IVSparse