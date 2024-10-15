/**
 * @file VCSC_Iterator_Methods.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Iterator Methods for VCSC Sparse Matrices
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    //* Constructors *//

    // Matrix Constructor
    template <typename T, typename indexT, bool columnMajor>
    inline SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::InnerIterator(IVSparse::SparseMatrix<T, indexT, 2, columnMajor>& matrix, uint32_t vec) {

        #ifdef IVSPARSE_DEBUG
        assert(vec < matrix.outerDim && vec >= 0 && "The vector index is out of bounds!");
        #endif

        // check if the vector is empty
        if (matrix.getNumUniqueVals(vec) == 0) {
            this->vals = nullptr;
            this->counts = nullptr;
            this->indices = nullptr;
            this->val = nullptr;
            this->index = 0;
            this->count = 0;

            countIndex = 0;
            indexSize = 0;

            return;
        }

        this->outer = vec;

        // set the pointers to the correct locations
        this->vals = matrix.getValues(vec);
        this->counts = matrix.getCounts(vec);
        this->indices = matrix.getIndices(vec);

        this->valsSize = matrix.getNumUniqueVals(vec);
        this->indexSize = matrix.getNumIndices(vec);

        // set the values of the iterator
        this->val = vals;
        this->index = indices[0];
        this->count = counts[0];
    }

    // Vector Constructor
    template <typename T, typename indexT, bool columnMajor>
    SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::InnerIterator(SparseMatrix<T, indexT, 2, columnMajor>::Vector& vector) {
        // check if the vector is empty
        if (vector.nonZeros() == 0) {
            this->vals = nullptr;
            this->counts = nullptr;
            this->indices = nullptr;
            this->val = nullptr;
            this->index = 0;
            this->count = 0;

            countIndex = 0;
            indexSize = 0;

            return;
        }

        this->outer = 0;

        // set the pointers to the correct locations
        this->vals = vector.getValues();
        this->counts = vector.getCounts();
        this->indices = vector.getIndices();

        this->valsSize = vector.uniqueVals();
        this->indexSize = vector.nonZeros();

        // set the values of the iterator
        this->val = vals;
        this->index = indices[0];
        this->count = counts[0];
    }

    //* Getters *//

    // Get the outer dimension
    template <typename T, typename indexT, bool columnMajor>
    indexT SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::outerDim() {
        return outer;
    }

    // Get the value
    template <typename T, typename indexT, bool columnMajor>
    T SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::value() {
        return *val;
    }

    // Get the index
    template <typename T, typename indexT, bool columnMajor>
    indexT SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::getIndex() {
        return index;
    }

    // Get a pointer to the value
    template <typename T, typename indexT, bool columnMajor>
    inline void SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::coeff(
        T newValue) {
        *val = newValue;
    }

    // Get the current row
    template <typename T, typename indexT, bool columnMajor>
    indexT SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::row() {
        if constexpr (!columnMajor) {
            return outer;
        }
        else {
            return index;
        }
    }

    // Get the current column
    template <typename T, typename indexT, bool columnMajor>
    indexT SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::col() {
        if constexpr (!columnMajor) {
            return index;
        }
        else {
            return outer;
        }
    }

    //* Operator Overloads *//

    // Dereference Operator
    template <typename T, typename indexT, bool columnMajor>
    T& SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::operator*() {
        return *val;
    }

    // Equality Operator
    template <typename T, typename indexT, bool columnMajor>
    bool SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::operator==(const InnerIterator& other) {
        return values == other.values;
    }

    // Inequality Operator
    template <typename T, typename indexT, bool columnMajor>
    bool SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::operator!=(const InnerIterator& other) {
        return values != other.values;
    }

    // Less Than Operator
    template <typename T, typename indexT, bool columnMajor>
    bool SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::operator<(const InnerIterator& other) {
        return values < other.values;
    }

    // Greater Than Operator
    template <typename T, typename indexT, bool columnMajor>
    bool SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::operator>(const InnerIterator& other) {
        return values > other.values;
    }

    // Increment Operator
    template <typename T, typename indexT, bool columnMajor>
    inline void SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator::operator++() {

        // decriment count
        count--;

        // if count is 0 then we need to move to the next value
        if (count == 0) {
            // Check if we are at the end of the values array
            countIndex++;

            if (countIndex >= indexSize) {
                return;
            }

            // Move to the next value
            val++;

            // Get the new count
            counts++;
            count = *counts;

            // Get the new index
            indices++;
            index = *indices;
        }
        else {
            // Move to the next index
            countIndex++;
            indices++;
            index = *indices;
        }

    }  // end operator++

}  // end namespace IVSparse