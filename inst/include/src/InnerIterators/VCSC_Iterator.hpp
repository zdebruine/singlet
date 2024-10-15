/**
 * @file VCSC_Iterator.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Inner Iterator for VCSC Declerations
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    /**
     * VCSC Inner Iterator Class \n \n
     * The VCSC Inner Iterator is a forward traversal iterator like the others in
     * the IVSparse library. It's very low overhead and is used to traverse over the
     * nonzeros of a single vector of a matrix or a vector on its own. The VCSC
     * Inner Iterator differs from the CSC Iterator in that it travereses a counts
     * vector in the VCSC Matrix instead of a outer pointers vector.
     */
    template <typename T, typename indexT, bool columnMajor>
    class SparseMatrix<T, indexT, 2, columnMajor>::InnerIterator {
        private:
        //* Private Class Variables *//

        indexT outer = 0;     // Outer dimension
        indexT index = 0;     // Current index
        indexT newIndex = 0;  // Next index
        T* val = nullptr;     // Current value

        T* vals = nullptr;          // Pointer to values
        indexT* counts = nullptr;   // Pointer to counts
        indexT* indices = nullptr;  // Pointer to indices

        indexT valsSize = 0;   // Number of unique values
        indexT indexSize = 0;  // Number of indices

        indexT count = 0;       // Current count
        indexT countIndex = 0;  // Current count of indices

        //* Private Class Methods *//

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
         * VCSC Matrix InnerIterator Constructor \n \n
         * The main constructor for the Inner Iterator. Given a matrix the iterator
         * will forward traverse over the given vector of the matrix. The traversal
         * is sorted by value in ascending order.
         */
        InnerIterator(SparseMatrix<T, indexT, 2, columnMajor>& mat, uint32_t col);

        /**
         * VCSC Vector InnerIterator Constructor \n \n
         * Same as the previous constructor but for a single standalone vector.
         * Can be used in the same way as the previous constructor.
         */
        InnerIterator(SparseMatrix<T, indexT, 2, columnMajor>::Vector& vec);

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
         */
        void coeff(T newValue);

        /**
         * @returns If the iterator is at the beginning of a new run.
         */

         ///@}

         //* Operator Overloads *//

         // Prefix increment operator
        void __attribute__((hot)) operator++();

        // Equality operator
        bool operator==(const InnerIterator& other);

        // Inequality operator
        bool operator!=(const InnerIterator& other);

        // Less than operator
        bool operator<(const InnerIterator& other);

        // Greater than operator
        bool operator>(const InnerIterator& other);

        // Boolean operator
        inline __attribute__((hot)) operator bool() { return countIndex < indexSize; }

        // Dereference operator
        T& operator*();

    };  // End of VCSC Inner Iterator Class

}  // namespace IVSparse