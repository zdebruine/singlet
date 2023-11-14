/**
 * @file IVCSC_SparseMatrix.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief IVCSC Sparse Matrix Class Declarations
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    /**
     * @tparam T The data type of the values in the matrix
     * @tparam indexT The data type of the indices in the matrix
     * @tparam compressionLevel The compression level used
     * @tparam columnMajor Whether the matrix is stored in column major format
     *
     * A class to represent a sparse matrix compressed in the Compressed Sparse
     * Fiber format (IVSparse). \n \n IVSparse Sparse Matrix is a read-only matrix
     * class optimized for sparse-dense computation in cases where values are highly
     * redundant. For such cases, sparse fiber storage can reduce memory footprint
     * by up to 50% compared to standard sparse compression. IVSparse also increases
     * the ability to further compress index arrays within each run. This default
     * templated version is for compression3 specifically. For compression level 1
     * and 2 there are template specializations.
     */
    template <typename T, typename indexT = uint64_t, uint8_t compressionLevel = 3, bool columnMajor = true>
    class SparseMatrix {
        private:
        //* The Matrix Data *//

        void** data = nullptr;         // The data of the matrix
        void** endPointers = nullptr;  // The pointers to the end of each column

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

        //* Private Methods *//

        // Compression Algorithm for going from CSC to VCSC or IVCSC
        template <typename T2, typename indexT2>
        void compressCSC(T2* vals, indexT2* innerIndices, indexT2* outerPointers);


        // Takes info about the value type and encodes it into a single uint32_t
        void encodeValueType();

        // Checks the value type matches the class template T
        void checkValueType();

        // Does checks on the class to ensure it is valid
        void userChecks();

        // Method to calcuate and set the byte size of the matrix in memory
        void calculateCompSize();

        // Scalar Multiplication
        inline IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>scalarMultiply(T scalar);

        // In Place Scalar Multiplication
        inline void inPlaceScalarMultiply(T scalar);

        // Matrix Vector Multiplication
        inline Eigen::Matrix<T, -1, 1> vectorMultiply(Eigen::Matrix<T, -1, 1>& vec);

        // Matrix Vector Multiplication 2 (with IVSparse Vector)
        inline Eigen::Matrix<T, -1, 1> vectorMultiply(typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector& vec);

        // Matrix Matrix Multiplication
        inline Eigen::Matrix<T, -1, -1> matrixMultiply(Eigen::Matrix<T, -1, -1> mat);

        // helper for ostream operator
        void print(std::ostream& stream);

        public:

        // Gets the number of rows in the matrix
        uint32_t rows() const { return numRows; }

        // Gets the number of columns in the matrix
        uint32_t cols() const { return numCols; }

        // Gets the inner dimension of the matrix
        uint32_t innerSize() const { return innerDim; }

        // Gets the outer dimension of the matrix
        uint32_t outerSize() const { return outerDim; }

        // Gets the number of non-zero elements in the matrix
        uint32_t nonZeros() const { return nnz; }

        // Gets the number of bytes needed to store the matrix
        size_t byteSize() const { return compSize; }


        //* Nested Subclasses *//

        // Vector Class for IVCSC Sparse Matrix
        class Vector;

        // Iterator Class for IVCSC Sparse Matrix
        class InnerIterator;

        //* Constructors *//
        /** @name Constructors
         */
         ///@{

         /**
          * Construct an empty IVSparse matrix \n \n
          * The matrix will have 0 rows and 0 columns and
          * will not be initialized with any values. All data
          * will be set to nullptr.
          *
          * @warning This constructor is not recommended for use as updating a IVSparse
          * matrix is not well supported.
          */
        SparseMatrix() {};

        // Private Helper Constructor for tranposing a IVSparse matrix
        SparseMatrix(std::unordered_map<T, std::vector<indexT>>* maps, uint32_t num_rows, uint32_t num_cols);


        /**
         * Empty Constructor \n \n
         * Takes in the number of rows and cols desired for an all zero matrix
         * of the specified size. All data will be set to nullptr.
         */
        SparseMatrix(uint32_t num_rows, uint32_t num_cols);

        /**
         * @param mat The Eigen Sparse Matrix to be compressed
         *
         * Eigen Sparse Matrix Constructor \n \n
         * This constructor takes an Eigen Sparse Matrix and compresses it into a
         * IVSparse matrix.
         */
        SparseMatrix(Eigen::SparseMatrix<T>& mat);

        /**
         * @param mat The Eigen Sparse Matrix to be compressed
         *
         * Eigen Sparse Matrix Constructor (Row Major) \n \n
         * Same as previous constructor but for Row Major Eigen Sparse Matrices.
         */
        SparseMatrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& mat);

        /**
         * @tparam compressionLevel2 The compression level of the IVSparse matrix to
         * convert
         * @param mat The IVSparse matrix to convert
         *
         * Convert a IVSparse matrix of a different compression level to this
         * compression level. \n \n This constructor takes in a IVSparse matrix of the
         * same storage order, value, and index type and converts it to a different
         * compresion level. This is useful for converting between compression levels
         * without having to go through the CSC format.
         */
        template <uint8_t compressionLevel2>
        SparseMatrix(IVSparse::SparseMatrix<T, indexT, compressionLevel2, columnMajor>& other);

        /**
         * @param other The IVSparse matrix to be copied
         *
         * Deep Copy Constructor \n \n
         * This constructor takes in a IVSparse matrix and creates a deep copy of it.
         */
        SparseMatrix(const IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>& other);

        /**
         * Raw CSC Constructor \n \n
         * This constructor takes in raw CSC storage format pointers and converts it
         * to a IVSparse matrix. One could also take this information and convert to
         * an Eigen Sparse Matrix and then to a IVSparse matrix.
         */
        template <typename T2, typename indexT2>
        SparseMatrix(T2* vals, indexT2* innerIndices, indexT2* outerPtr, uint32_t num_rows, uint32_t num_cols, uint32_t nnz);

        /**
         * COO Tuples Constructor \n \n
         * This constructor takes in a list of tuples in COO format which can be
         * unsorted but without duplicates. The tuples are sorted and then converted
         * to a IVSparse matrix.
         *
         * @note COO is (row, col, value) format.
         *
         * @warning This constructor does not allow for duplicates but will sort the
         * tuples.
         */
        template <typename T2, typename indexT2>
        SparseMatrix(std::vector<std::tuple<indexT2, indexT2, T2>>& entries, uint64_t num_rows, uint32_t num_cols, uint32_t nnz);

        /**
         * @param vec The vector to construct the matrix from
         *
         * IVSparse Vector Constructor \n \n
         * This constructor takes in a single IVSparse vector and creates a one
         * column/row IVSparse matrix.
         */
        SparseMatrix(typename IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector& vec);

        /**
         * @param vecs The vector of IVSparse vectors to construct from.
         *
         * Vector of IVSparse Vectors Constructor \n \n
         * This constructor takes in an vector of IVSparse vectors and creates a
         * IVSparse matrix from them.
         */
        SparseMatrix(std::vector<typename IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector>& vecs);

        /**
         * @param filename The filepath of the matrix to be read in
         *
         * File Constructor \n \n
         * Given a filepath to a IVSparse matrix written to file this constructor will
         * read in the matrix and construct it.
         */
        SparseMatrix(const char* filename);

        /**
         * @brief Destroy the Sparse Matrix object
         */
        ~SparseMatrix();

        ///@}

        //* Getters *//
        /**
         * @name Getters
         */
         ///@{

         /**
          * @returns T The value at the specified row and column. Returns 0 if the
          * value is not found.
          *
          * Get the value at the specified row and column
          *
          * @note Users cannot update individual values in a IVSparse matrix.
          *
          * @warning This method is not efficient and should not be used in performance
          * critical code.
          */
        T coeff(uint32_t row, uint32_t col);

        /**
         * @returns true If the matrix is stored in column major format
         * @returns false If the matrix is stored in row major format
         *
         * See the storage order of the IVSparse matrix.
         */
        bool isColumnMajor() const;

        /**
         * @param vec The vector to get the pointer to
         * @returns void* The pointer to the vector
         *
         * Get a pointer to a vector in the IVSparse matrix such as the first column.
         *
         * @note Can only get vectors in the storage order of the matrix.
         */
        void* vectorPointer(uint32_t vec);

        /**
         * @param vec The vector to get a copy of
         * @returns Vector The vector copy returned
         *
         * Get a copy of a IVSparse vector from the IVSparse matrix such as the first
         * column.
         *
         * @note Can only get vectors in the storage order of the matrix.
         */
        typename IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector getVector(uint32_t vec);

        /**
         * @param vec The vector to get the size of
         * @returns size_t The size of the vector in bytes
         *
         * Get the size of a vector in the IVSparse matrix in bytes.
         *
         * @note Can only get vectors in the storage order of the matrix.
         */
        size_t getVectorSize(uint32_t vec) const;

        ///@}

        //* Calculations *//
        /**
         * @name Calculations
         */
         ///@{

         /**
          * @returns A vector of the sum of each vector along the outer dimension.
          */
        inline std::vector<T> outerSum();

        /**
         * @returns A vector of the sum of each vector along the inner dimension.
         */
        inline std::vector<T> innerSum();

        /**
         * @returns A vector of the maximum value in each column.
         */
        inline std::vector<T> maxColCoeff();

        /**
         * @returns A vector of the maximum value in each row.
         */
        inline std::vector<T> maxRowCoeff();

        /**
         * @returns A vector of the minimum value in each column.
         */
        inline std::vector<T> minColCoeff();

        /**
         * @returns A vector of the minimum value in each row.
         */
        inline std::vector<T> minRowCoeff();

        /**
         * @returns The trace of the matrix.
         *
         * @note Only works for square matrices.
         */
        inline T trace();

        /**
         * @returns The sum of all the values in the matrix.
         */
        inline T sum();

        /**
         * @returns The frobenius norm of the matrix.
         */
        inline double norm();

        /**
         * @returns Returns the length of the specified vector.
         */
        inline double vectorLength(uint32_t vec);

        ///@}

        //* Utility Methods *//
        /**
         * @name Utility Methods
         */
         ///@{

         /**
          * @param filename The filename of the matrix to write to
          *
          * This method writes the IVSparse matrix to a file in binary format.
          * This can then be read in later using the file constructor.
          * Currently .ivsparse is the perfered file extension.
          *
          * @note Useful to split a matrix up and then write each part separately.
          */
        void write(const char* filename);

        /**
         * Prints "IVSparse Matrix:" followed by the dense representation of the
         * matrix to the console.
         *
         * @note Useful for debugging but only goes up to 100 of either dimension.
         */
        void print();

        /**
         * @returns The current matrix as uncompressed to CSC format.
         */
        IVSparse::SparseMatrix<T, indexT, 1, columnMajor> toCSC();

        /**
         * @returns The current matrix as a VCSC Matrix.
         */
        IVSparse::SparseMatrix<T, indexT, 2, columnMajor> toVCSC();

        /**
         * @returns An Eigen Sparse Matrix constructed from the IVSparse matrix data.
         */
        Eigen::SparseMatrix<T, columnMajor ? Eigen::ColMajor : Eigen::RowMajor> toEigen();

        ///@}

        //* Matrix Manipulation Methods *//
        /**
         * @name Matrix Manipulation Methods
         */
         ///@{

         /**
          * @returns A transposed version of the IVSparse matrix.
          *
          * @warning This method is not very efficient for VCSC and IVCSC matrices.
          */
        IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor> transpose();

        /**
         * Transposes the matrix in place instead of returning a new matrix.
         *
         * @warning This method is not very efficient for VCSC and IVCSC matrices.
         */
        void inPlaceTranspose();

        /**
         * @param mat The matrix to append to the matrix in the correct storage order.
         *
         * Appends an IVSparse matrix to the current matrix in the storage order of the
         * matrix.
         */
        void append(SparseMatrix<T, indexT, compressionLevel, columnMajor>& mat);

        /**
         * @param mat The matrix to append to the matrix in the correct storage order.
         *
         * Appends an Eigen::SparseMatrix to the current matrix in the storage order of the
         * matrix. This converts the Eigen::SparseMatrix to an IVSparse matrix.
         */

        inline void append(Eigen::SparseMatrix<T>& mat);


        /**
         * @returns A matrix that represent a slice of the
         * IVSparse matrix.
         */
        IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor> slice(uint32_t start, uint32_t end);

        ///@}

        //* Operator Overloads *//

        friend std::ostream& operator<< (std::ostream& stream, IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>& mat) {
            mat.print(stream);
            return stream;
        }

        // Assignment Operator
        IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>& operator=(const IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>& other);

        // Equality Operator
        bool operator==(const SparseMatrix<T, indexT, compressionLevel, columnMajor>& other) const;

        // Inequality Operator
        bool operator!=(const SparseMatrix<T, indexT, compressionLevel, columnMajor>& other);

        // Coefficient Access Operator
        T operator()(uint32_t row, uint32_t col);

        // Vector Access Operator
        typename IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector operator[](uint32_t vec);

        // Scalar Multiplication
        IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor> operator*(T scalar);

        // In Place Scalar Multiplication
        void operator*=(T scalar);

        // Matrix Vector Multiplication
        Eigen::Matrix<T, -1, 1> operator*(Eigen::Matrix<T, -1, 1>& vec);

        // Matrix Vector Multiplication 2 (with IVSparse Vector)
        Eigen::Matrix<T, -1, 1> operator*(typename SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector& vec);

        // Matrix Matrix Multiplication
        Eigen::Matrix<T, -1, -1> operator*(Eigen::Matrix<T, -1, -1>& mat);

    };  // End of SparseMatrix Class

}  // namespace IVSparse