/**
 * @file CSC_Iterator_Methods.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Iterator Methods for CSC Sparse Matrices
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    //* Constructors *//

    // CSC Matrix Constructor
    template <typename T, typename indexT, bool columnMajor>
    inline SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::InnerIterator(
        SparseMatrix<T, indexT, 1, columnMajor>& mat, uint32_t vec) {

        this->outer = vec;

        // check if the vector is empty
        if (mat.getOuterPointers()[vec] == mat.getOuterPointers()[vec + 1]) {
            vals = nullptr;
            indices = nullptr;
            endPtr = nullptr;
            return;
        }

        // set the pointers to the correct locations
        vals = &mat.vals[mat.outerPtr[vec]];
        indices = &mat.innerIdx[mat.outerPtr[vec]];
        endPtr = &mat.innerIdx[mat.outerPtr[vec + 1]];

        // set the values of the iterator
        val = vals;
        index = indices[0];
    }

    // CSC Vector Constructor
    template <typename T, typename indexT, bool columnMajor>
    SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::InnerIterator(
        SparseMatrix<T, indexT, 1, columnMajor>::Vector& vec) {

        this->outer = 0;

        // set the pointers to the correct locations
        vals = vec.values();
        indices = vec.indexPtr();
        endPtr = vec.indexPtr() + vec.nonZeros();

        // set the values of the iterator
        val = vals;
        index = indices[0];
    }

    //* Overloaded Operators *//

    // Increment Operator
    template <typename T, typename indexT, bool columnMajor>
    inline void SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::operator++() {
        vals++;
        indices++;

        // check if the iterator is at the end of the vector
        if (indices == endPtr) {
            return;
        }

        // set the values of the iterator
        val = vals;
        index = *indices;
    }

    // Equality Operator
    template <typename T, typename indexT, bool columnMajor>
    bool SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::operator==(const InnerIterator& other) {
        return (vals == other.vals && indices == other.index);
    }

    // Inequality Operator
    template <typename T, typename indexT, bool columnMajor>
    bool SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::operator!=(const InnerIterator& other) {
        return (vals != other.vals || indices != other.index);
    }

    // Less Than Operator
    template <typename T, typename indexT, bool columnMajor>
    bool SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::operator<(const InnerIterator& other) {
        return (vals < other.vals && indices < other.index);
    }

    // Greater Than Operator
    template <typename T, typename indexT, bool columnMajor>
    bool SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::operator>(const InnerIterator& other) {
        return (vals > other.vals && indices > other.index);
    }

    // Dereference Operator
    template <typename T, typename indexT, bool columnMajor>
    T& SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::operator*() {
        return val;
    }

    //* Getters & Setters *//

    // Get the current index of the iterator
    template <typename T, typename indexT, bool columnMajor>
    indexT SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::getIndex() {
        return index;
    }

    // Get the current outer dimension of the iterator
    template <typename T, typename indexT, bool columnMajor>
    indexT SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::outerDim() {
        return outer;
    }

    // Get the current row of the iterator
    template <typename T, typename indexT, bool columnMajor>
    indexT SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::row() {
        if (columnMajor) {
            return index;
        }
        else {
            return outer;
        }
    }

    // Get the current column of the iterator
    template <typename T, typename indexT, bool columnMajor>
    indexT SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::col() {
        if (columnMajor) {
            return outer;
        }
        else {
            return index;
        }
    }

    // Get the current value of the iterator
    template <typename T, typename indexT, bool columnMajor>
    T SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::value() {
        return *val;
    }

    // coefficent access method
    template <typename T, typename indexT, bool columnMajor>
    void SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator::coeff(T value) {
        *val = value;
    }

}  // namespace IVSparse