/**
 * @file CSC_Iterator.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Inner Iterator for CSC Declerations
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    /**
     * CSC Inner Iterator Class \n \n
     * The CSC Inner Iterator is a forward traversal iterator like the others in the
     * IVSparse library. It's very low overhead and is used to traverse over the
     * nonzeros of a single vector of a matrix or a vector on its own.
     */
    template <typename T, typename indexT, bool columnMajor>
    class SparseMatrix<T, indexT, 1, columnMajor>::InnerIterator {
        private:
        //* Private Class Variables *//

        T* val;        // Current value
        indexT index;  // Current index
        indexT outer;  // Outer dimension

        T* vals;
        indexT* indices;
        indexT* endPtr;

        public:
        //* Constructors & Destructor *//
        /** @name Constructors
         */
         ///@{

         /**
          * Default Iterator Constructor \n \n
          * Creates an empty iterator that can't be used on its own.
          */
        InnerIterator() {};

        /**
         * CSC Matrix InnerIterator Constructor \n \n
         * The main constructor for the Inner Iterator. Given a matrix the iterator
         * will forward traverse over the given vector of the matrix. The traversal
         * is sorted by index.
         */
        InnerIterator(SparseMatrix<T, indexT, 1, columnMajor>& mat, uint32_t vec);

        /**
         * CSC Vector InnerIterator Constructor \n \n
         * Same as the previous constructor but for a single standalone vector.
         * Can be used in the same way as the previous constructor.
         */
        InnerIterator(SparseMatrix<T, indexT, 1, columnMajor>::Vector& vec);

        ///@}

        //* Getters *//
        /** @name Getters
         */
         ///@{

         /**
          * @returns The current index of the iterator.
          */
        indexT getIndex();

        /**
         * @returns The current outer dimension of the iterator.
         */
        indexT outerDim();

        /**
         * @returns The current row of the iterator.
         */
        indexT row();

        /**
         * @returns The current column of the iterator.
         */
        indexT col();

        /**
         * @returns The current value of the iterator.
         */
        T value();

        /**
         * Changes the value where the iterator is pointing.
         *
         * @note This is the only way to update elements in the IVSparse format.
         *
         * @warning This method may break things if used without care, IVSparse is not
         * meant to update values.
         */
        void coeff(T newValue);

        ///@}

        //* Operator Overloads *//

        // Prefix Increment Operator
        void __attribute__((hot)) operator++();

        // Equality Operator
        bool operator==(const InnerIterator& other);

        // Inequality Operator
        bool operator!=(const InnerIterator& other);

        // Less Than Operator
        bool operator<(const InnerIterator& other);

        // Greater Than Operator
        bool operator>(const InnerIterator& other);

        // Dereference Operator
        T& operator*();

        // Bool Operator
        inline __attribute__((hot)) operator bool() { return indices < endPtr; };

    };  // class InnerIterator

}  // namespace IVSparse