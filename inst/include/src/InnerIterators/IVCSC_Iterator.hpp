/**
 * @file IVCSC_Iterator.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Inner Iterator for IVCSC Declerations
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    /**
     * @tparam T The type of the values in the matrix
     * @tparam indexT The type of the indices in the matrix
     * @tparam compressionLevel The level of compression used in the matrix
     * @tparam columnMajor Whether the matrix is column major or not
     *
     * IVCSC Inner Iterator Class \n \n
     * The IVCSC Inner Iterator is a forward traversal iterator like the others in
     * the IVSparse library. The IVCSC Iterator is slower than the others due to
     * needing to decode compressed data.
     */
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    class SparseMatrix<T, indexT, compressionLevel, columnMajor>::InnerIterator {
        private:
        //* Private Class Variables *//

        indexT outer;      // Outer dimension
        indexT index;      // Current index
        T* val = nullptr;  // Current value

        indexT newIndex = 0;  // Next index

        uint8_t indexWidth = 1;  // Width of the current run

        void* data;    // Pointer to the current data
        void* endPtr;  // Pointer to the end of the data

        bool firstIndex = true;  // Is this the first index of the vector

        //* Private Class Methods *//

        // Decodes the index from the data pointer
        void __attribute__((hot)) decodeIndex();

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
         * IVCSC Matrix InnerIterator Constructor \n \n
         * The main constructor for the Inner Iterator. Given a matrix the iterator
         * will forward traverse over the given vector of the matrix. The traversal
         * is sorted by value in ascending order.
         */
        InnerIterator(SparseMatrix<T, indexT, compressionLevel, columnMajor>& mat,
                      uint32_t col);

        /**
         * IVCSC Vector InnerIterator Constructor \n \n
         * Same as the previous constructor but for a single standalone vector.
         * Can be used in the same way as the previous constructor.
         */
        InnerIterator(
            SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector& vec);

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
        bool isNewRun();

        ///@}

        //* Operator Overloads *//

        // Increment Operator
        void __attribute__((hot)) operator++();

        // Equality Operators
        bool operator==(const InnerIterator& other);

        // Inequality Operators
        bool operator!=(const InnerIterator& other);

        // Less Than Operator
        bool operator<(const InnerIterator& other);

        // Greater Than Operator
        bool operator>(const InnerIterator& other);

        // Bool Operator
        inline __attribute__((hot)) operator bool() {
            return ((char*)endPtr - indexWidth > data);
        }

        // Dereference Operator
        T& operator*();

    };  // End of InnerIterator Class

}  // namespace IVSparse
