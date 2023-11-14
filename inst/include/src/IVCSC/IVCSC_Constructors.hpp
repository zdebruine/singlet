/**
 * @file IVCSC_Constructors.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Constructors for IVCSC Sparse Matrices
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

    // Destructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::~SparseMatrix() {
        // delete the meta data
        if (metadata != nullptr) {
            delete[] metadata;
        }

        if (data != nullptr) {
            // free the data
            for (size_t i = 0; i < outerDim; i++) {
                if (data[i] != nullptr) {
                    free(data[i]);
                }
            }
            free(data);
        }

        if (endPointers != nullptr) {
            free(endPointers);
        }
    }

    // Row and Column Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(uint32_t num_rows, uint32_t num_cols) {

        #ifdef IVSPARSE_DEBUG
        // check that the number of rows and columns is greater than 0
        assert(num_rows > 0 && num_cols > 0 &&
               "The number of rows and columns must be greater than 0");
        #endif

        // set class variables
        if constexpr (columnMajor) {
            innerDim = num_cols;
            outerDim = num_rows;
        }
        else {
            innerDim = num_rows;
            outerDim = num_cols;
        }
        numRows = num_rows;
        numCols = num_cols;
        encodeValueType();
        index_t = sizeof(indexT);
        nnz = 0;

        // allocate memory for the data
        try {
            data = (void**)malloc(outerDim * sizeof(void*));
            endPointers = (void**)malloc(outerDim * sizeof(void*));
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << '\n';
        }

        // set all data and endpointers to the nullptr
        for (size_t i = 0; i < outerDim; i++) {
            data[i] = nullptr;
            endPointers[i] = nullptr;
        }

        // set the metadata
        metadata = new uint32_t[NUM_META_DATA];
        metadata[0] = compressionLevel;
        metadata[1] = innerDim;
        metadata[2] = outerDim;
        metadata[3] = nnz;
        metadata[4] = val_t;
        metadata[5] = index_t;

        // calculate the compression size
        calculateCompSize();

        // run the user checks
        #ifdef IVSPARSE_DEBUG
        userChecks();
        #endif
    }

    // Eigen Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(Eigen::SparseMatrix<T>& mat) {

        // check if the matrix is empty
        if (mat.nonZeros() == 0) {
            val_t = 0;
            index_t = 0;

            data = nullptr;
            endPointers = nullptr;
            metadata = nullptr;
            return;
        }

        // make sure the matrix is compressed
        mat.makeCompressed();

        // get the number of rows and columns
        numRows = mat.rows();
        numCols = mat.cols();

        outerDim = columnMajor ? numCols : numRows;
        innerDim = columnMajor ? numRows : numCols;

        // get the number of non-zero elements
        nnz = mat.nonZeros();

        // call the compression function
        compressCSC(mat.valuePtr(), mat.innerIndexPtr(), mat.outerIndexPtr());
    }

    // Eigen Row Major Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& mat) {

        // check if the matrix is empty
        if (mat.nonZeros() == 0) {
            val_t = 0;
            index_t = 0;

            data = nullptr;
            endPointers = nullptr;
            metadata = nullptr;
            return;
        }

        // make sure the matrix is compressed
        mat.makeCompressed();

        // get the number of rows and columns
        numRows = mat.rows();
        numCols = mat.cols();

        outerDim = numRows;
        innerDim = numCols;

        // get the number of non-zero elements
        nnz = mat.nonZeros();

        // call the compression function
        compressCSC(mat.valuePtr(), mat.innerIndexPtr(), mat.outerIndexPtr());
    }

    // Deep Copy Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(const IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>& other) {

        *this = other;
    }

    // Conversion Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    template <uint8_t otherCompressionLevel>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(IVSparse::SparseMatrix<T, indexT, otherCompressionLevel, columnMajor>& other) {

        // if already the right compression level
        if constexpr (otherCompressionLevel == compressionLevel) {
            *this = other;
            return;
        }

        // make a temporary matrix of the correct compression level
        IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor> temp;

        // convert other to the right compression level
        if constexpr (otherCompressionLevel == 1) {
            temp = other.toIVCSC();
        }
        else if constexpr (otherCompressionLevel == 2) {
            temp = other.toIVCSC();
        }

        // other should be the same compression level as this now
        *this = temp;

        // run the user checks
        #ifdef IVSPARSE_DEBUG
        userChecks();
        #endif
    }

    // Raw CSC Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    template <typename T2, typename indexT2>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(T2* vals, indexT2* innerIndices, indexT2* outerPtr, uint32_t num_rows, uint32_t num_cols, uint32_t nnz) {

        #ifdef IVSPARSE_DEBUG
        assert(num_rows > 0 && num_cols > 0 &&
               "The number of rows and columns must be greater than 0");
        assert(nnz > 0 && "The number of non-zero elements must be greater than 0");
        #endif

        // see if the matrix is empty
        if (nnz == 0) [[unlikely]] {
            val_t = 0;
            index_t = 0;
            data = nullptr;
            endPointers = nullptr;
            metadata = nullptr;
            return;
            }

            // set class variables
                if (columnMajor) {
                    innerDim = num_rows;
                    outerDim = num_cols;
                }
                else {
                    innerDim = num_cols;
                    outerDim = num_rows;
                }
            numRows = num_rows;
            numCols = num_cols;
            this->nnz = nnz;

            // call the compression function
            compressCSC(vals, innerIndices, outerPtr);
    }

    // COO Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    template <typename T2, typename indexT2>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(std::vector<std::tuple<indexT2, indexT2, T2>>& entries, uint64_t num_rows, uint32_t num_cols, uint32_t nnz) {


        #ifdef IVSPARSE_DEBUG
        assert(num_rows > 0 && num_cols > 0 &&
               "The number of rows and columns must be greater than 0");
        assert(nnz > 0 && "The number of non-zero elements must be greater than 0");
        #endif

        // see if the matrix is empty
        if (nnz == 0) [[unlikely]] {
            *this = SparseMatrix<T, indexT, compressionLevel, columnMajor>(num_rows, num_cols);
            return;
            }

            // set class variables
                if (columnMajor) {
                    innerDim = num_rows;
                    outerDim = num_cols;
                }
                else {
                    innerDim = num_cols;
                    outerDim = num_rows;
                }

            numRows = num_rows;
            numCols = num_cols;
            this->nnz = nnz;
            encodeValueType();
            index_t = sizeof(indexT);

            metadata = new uint32_t[NUM_META_DATA];
            metadata[0] = compressionLevel;
            metadata[1] = innerDim;
            metadata[2] = outerDim;
            metadata[3] = nnz;
            metadata[4] = val_t;
            metadata[5] = index_t;

            // allocate memory for the data
            try {
                data = (void**)malloc(outerDim * sizeof(void*));
                endPointers = (void**)malloc(outerDim * sizeof(void*));
            }
            catch (std::bad_alloc& e) {
                std::cerr << "Error: Could not allocate memory for IVSparse matrix"
                    << std::endl;
                exit(1);
            }

            // set all data and endpointers to the nullptr
            for (size_t i = 0; i < outerDim; i++) {
                data[i] = nullptr;
                endPointers[i] = nullptr;
            }

            // sort the tuples by first by column then by row
            std::sort(entries.begin(), entries.end(),
                      [](const std::tuple<indexT2, indexT2, T2>& a,
                         const std::tuple<indexT2, indexT2, T2>& b) {
                             if (std::get<1>(a) == std::get<1>(b)) {
                                 return std::get<0>(a) < std::get<0>(b);
                             }
                             else {
                                 return std::get<1>(a) < std::get<1>(b);
                             }
                      });

            std::map<T2, std::vector<indexT2>> maps[outerDim];

            // loop through the tuples
            for (size_t i = 0; i < nnz; i++) {
                // get the column
                indexT2 row = std::get<0>(entries[i]);
                indexT2 col = std::get<1>(entries[i]);
                T2 val = std::get<2>(entries[i]);

                // check if the value is already in the map
                if (maps[col].find(val) != maps[col].end()) {
                    // value found positive delta encode it
                    maps[col][val].push_back(row - maps[col][val][1]);

                    // update the last index
                    maps[col][val][1] = row;

                    // update the maximum delta
                    if (maps[col][val][maps[col][val].size() - 1] > maps[col][val][0])
                        maps[col][val][0] = maps[col][val][maps[col][val].size() - 1];
                }
                else {
                    // value not found
                    maps[col][val] = std::vector<indexT2>{ row };

                    // add maximum delta and last index placeholders
                    maps[col][val].push_back(row);
                    maps[col][val].push_back(row);
                }
            }  // end of loop through tuples

            // loop through the maps
            #ifdef IVSPARSE_HAS_OPENMP
            #pragma omp parallel for
            #endif
            for (uint32_t i = 0; i < outerDim; i++) {
                size_t outerByteSize = 0;

                for (auto& pair : maps[i]) {
                    // change first value to be byte width of the maximum delta
                    pair.second[0] = byteWidth(pair.second[0]);

                    // add the size of the run to the size of the column
                    //* value + index width + indices * index width + delimiter (index width)
                    outerByteSize += sizeof(T) + 1 +
                        (pair.second[0] * (pair.second.size() - 2)) +
                        pair.second[0];
                }

                // if column is empty set the data and endpointer to nullptr
                if (outerByteSize == 0) {
                    data[i] = nullptr;
                    endPointers[i] = nullptr;
                    continue;
                }

                // allocate space for the column
                try {
                    data[i] = malloc(outerByteSize);
                }
                catch (std::bad_alloc& e) {
                    std::cout << "Error: " << e.what() << std::endl;
                    exit(1);
                }

                // get a help pointer for moving through raw memory
                void* helpPtr = data[i];

                // loop through the dictionary and write to memory
                for (auto& pair : maps[i]) {
                    // Write the value to memory
                    *(T*)helpPtr = (T)pair.first;
                    helpPtr = (T*)helpPtr + 1;

                    // also write the index width
                    *(uint8_t*)helpPtr = (uint8_t)pair.second[0];
                    helpPtr = (uint8_t*)helpPtr + 1;

                    // loop through the indices and write them to memory
                    for (size_t k = 2; k < pair.second.size(); k++) {
                        // if compression level 3 skip the first two indices and cast the index


                        // create a type of the correct width
                        switch (pair.second[0]) {
                        case 1:
                            *(uint8_t*)helpPtr = (uint8_t)pair.second[k];
                            helpPtr = (uint8_t*)helpPtr + 1;
                            break;
                        case 2:
                            *(uint16_t*)helpPtr = (uint16_t)pair.second[k];
                            helpPtr = (uint16_t*)helpPtr + 1;
                            break;
                        case 3:
                            *(uint32_t*)helpPtr = (uint32_t)pair.second[k] & 0xFFFFFF;
                            helpPtr = (uint8_t*)helpPtr + 3;
                            break;
                        case 4:
                            *(uint32_t*)helpPtr = (uint32_t)pair.second[k];
                            helpPtr = (uint32_t*)helpPtr + 1;
                            break;
                        case 5:
                            *(uint64_t*)helpPtr = (uint64_t)pair.second[k] & 0xFFFFFFFFFF;
                            helpPtr = (uint8_t*)helpPtr + 5;
                            break;
                        case 6:
                            *(uint64_t*)helpPtr = (uint64_t)pair.second[k] & 0xFFFFFFFFFFFF;
                            helpPtr = (uint8_t*)helpPtr + 6;
                            break;
                        case 7:
                            *(uint64_t*)helpPtr = (uint64_t)pair.second[k] & 0xFFFFFFFFFFFFFF;
                            helpPtr = (uint8_t*)helpPtr + 7;
                            break;
                        case 8:
                            *(uint64_t*)helpPtr = (uint64_t)pair.second[k];
                            helpPtr = (uint64_t*)helpPtr + 1;
                            break;
                        }

                    }  // End of index loop

                    // write a delimiter of the correct width
                    memset(helpPtr, 0, pair.second[0]);
                    helpPtr = (uint8_t*)helpPtr + pair.second[0];

                    // Set a pointer to the end of the data
                    endPointers[i] = helpPtr;

                }  // End of dictionary loop
            }

            // calculate the compression size
            calculateCompSize();
    }

    // IVSparse Vector Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(typename IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector& vec) {


        // if it is the matrix is empty and we need to construct it
        if (columnMajor) {
            numRows = vec.getLength();
            numCols = 1;
            innerDim = numRows;
            outerDim = numCols;
        }
        else {
            numRows = 1;
            numCols = vec.getLength();
            innerDim = numCols;
            outerDim = numRows;
        }

        nnz = vec.nonZeros();

        encodeValueType();
        index_t = sizeof(indexT);

        // malloc the data
        try {
            data = (void**)malloc(sizeof(void*));
            endPointers = (void**)malloc(sizeof(void*));
        }
        catch (std::bad_alloc& e) {
            throw std::bad_alloc();
        }

        // malloc the vector
        try {
            // if vector is empty set the data to null
            if (vec.begin() == vec.end()) {
                data[0] = nullptr;
                endPointers[0] = nullptr;
            }
            else {
                data[0] = malloc(vec.byteSize());
                endPointers[0] = (char*)data[0] + vec.byteSize();
            }
        }
        catch (std::bad_alloc& e) {
            throw std::bad_alloc();
        }

        // copy the vector into the matrix if the vector is not empty
        if (vec.begin() != vec.end()) memcpy(data[0], vec.begin(), vec.byteSize());

        // set the metadata
        metadata = new uint32_t[NUM_META_DATA];
        metadata[0] = compressionLevel;
        metadata[1] = innerDim;
        metadata[2] = outerDim;
        metadata[3] = nnz;
        metadata[4] = val_t;
        metadata[5] = index_t;

        // run the user checks and calculate the compression size
        calculateCompSize();
        #ifdef IVSPARSE_DEBUG
        userChecks();
        #endif
    }

    // Array of Vectors Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(std::vector<typename IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor>::Vector>& vecs) {

        IVSparse::SparseMatrix<T, indexT, compressionLevel, columnMajor> temp(vecs[0]);

        for (size_t i = 1; i < vecs.size(); i++) {
            temp.append(vecs[i]);
        }

        *this = temp;

        calculateCompSize();

        // run the user checks
        #ifdef IVSPARSE_DEBUG
        userChecks();
        #endif
    }

    // File Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(const char* filename) {

        FILE* fp = fopen(filename, "rb");

        #ifdef IVSPARSE_DEBUG
        if (fp == nullptr) {
            throw std::runtime_error("Error: Could not open file");
        }
        #endif

        // read the metadata
        metadata = new uint32_t[NUM_META_DATA];
        if(fread(metadata, sizeof(uint32_t), NUM_META_DATA, fp) != 1) [[unlikely]] {
            throw std::runtime_error("Error: Could not read file");
        }

        // set the matrix info
        innerDim = metadata[1];
        outerDim = metadata[2];
        nnz = metadata[3];
        val_t = metadata[4];
        index_t = metadata[5];

        numRows = columnMajor ? innerDim : outerDim;
        numCols = columnMajor ? outerDim : innerDim;

        #ifdef IVSPARSE_DEBUG
        // if the compression level of the file is different than the compression
        // level of the class
        if (metadata[0] != compressionLevel) {
            // throw an error
            throw std::runtime_error(
                "Error: Compression level of file does not match compression level of "
                "class");
        }
        #endif

        // allocate the memory
        try {
            data = (void**)malloc(outerDim * sizeof(void*));
            endPointers = (void**)malloc(outerDim * sizeof(void*));
        }
        catch (std::bad_alloc& e) {
            std::cerr << "Error: Could not allocate memory for IVSparse matrix"
                << std::endl;
            exit(1);
        }

        // get the vector sizes and allocate memory
        for (size_t i = 0; i < outerDim; i++) {
            // get the size of the column
            uint64_t size;

            if(fread(&size, sizeof(uint64_t), 1, fp) != 1) [[unlikely]] {
                throw std::runtime_error("Error: Could not read file");
            }

            // if the size is 0, set the data and endpointer to nullptr
            if (size == 0) {
                data[i] = nullptr;
                endPointers[i] = nullptr;
                continue;
            }

            // malloc the column
            try {
                data[i] = malloc(size);
                endPointers[i] = (char*)data[i] + size;
            }
            catch (std::bad_alloc& e) {
                throw std::bad_alloc();
            }
        }

        // read the data
        for (size_t i = 0; i < outerDim; i++) {
            if(fread(data[i], 1, (uint8_t*)endPointers[i] - (uint8_t*)data[i], fp) != 1) [[unlikely]] {
                throw std::runtime_error("Error: Could not read file");
            }
        }

        // close the file
        fclose(fp);

        // calculate the compresssion size
        calculateCompSize();

        // run the user checks
        #ifdef IVSPARSE_DEBUG
        userChecks();
        #endif

    }  // end of file constructor

    //* Private Constructors *//

    // Private Tranpose Constructor
    template <typename T, typename indexT, uint8_t compressionLevel, bool columnMajor>
    SparseMatrix<T, indexT, compressionLevel, columnMajor>::SparseMatrix(std::unordered_map<T, std::vector<indexT>>* maps, uint32_t num_rows, uint32_t num_cols) {

        // set class variables
        if constexpr (columnMajor) {
            innerDim = num_cols;
            outerDim = num_rows;
        }
        else {
            innerDim = num_rows;
            outerDim = num_cols;
        }

        numRows = num_cols;
        numCols = num_rows;
        encodeValueType();
        index_t = sizeof(indexT);

        // allocate memory for the data
        try {
            data = (void**)malloc(outerDim * sizeof(void*));
            endPointers = (void**)malloc(outerDim * sizeof(void*));
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << '\n';
        }

        // set all data and endpointers to the nullptr
        for (size_t i = 0; i < outerDim; i++) {
            data[i] = nullptr;
            endPointers[i] = nullptr;
        }

        //* logic here

        // loop through the array
        #ifdef IVSPARSE_HAS_OPENMP
        #pragma omp parallel for
        #endif
        for (size_t i = 0; i < outerDim; i++) {
            // check if the column is empty
            if (maps[i].empty()) [[unlikely]] {
                data[i] = nullptr;
                endPointers[i] = nullptr;
                continue;
            }

            size_t byteSize = 0;

            // loop through the vectors of the map
            for (auto& val : maps[i]) {
                // add the size of the vector to the byteSize
                byteSize += sizeof(T) + 1 +
                    (val.second[val.second.size() - 1] * (val.second.size() - 1) +
                     val.second[val.second.size() - 1]);
            }

            // allocate memory for the vector
            try {
                data[i] = malloc(byteSize);
            }
            catch (const std::exception& e) {
                std::cerr << e.what() << '\n';
            }

            // set the end pointer
            endPointers[i] = (char*)data[i] + byteSize;

            // compressCSC the column
            void* helpPtr = data[i];

            for (auto& val : maps[i]) {
                nnz += val.second.size() - 1;

                // set the value
                *(T*)helpPtr = val.first;
                helpPtr = (char*)helpPtr + sizeof(T);
                *(uint8_t*)helpPtr = (uint8_t)val.second[val.second.size() - 1];
                helpPtr = (uint8_t*)helpPtr + 1;

                // write the indices
                for (size_t k = 0; k < val.second.size(); k++) {
                    if (k == val.second.size() - 1) break;

                    switch (val.second[val.second.size() - 1]) {
                    case 1:
                        *(uint8_t*)helpPtr = (uint8_t)val.second[k];
                        helpPtr = (uint8_t*)helpPtr + 1;
                        break;
                    case 2:
                        *(uint16_t*)helpPtr = (uint16_t)val.second[k];
                        helpPtr = (uint16_t*)helpPtr + 1;
                        break;
                    case 3:
                        *(uint32_t*)helpPtr = (uint32_t)val.second[k] & 0xFFFFFF;
                        helpPtr = (uint8_t*)helpPtr + 3;
                        break;
                    case 4:
                        *(uint32_t*)helpPtr = (uint32_t)val.second[k];
                        helpPtr = (uint32_t*)helpPtr + 1;
                        break;
                    case 5:
                        *(uint64_t*)helpPtr = (uint64_t)val.second[k] & 0xFFFFFFFFFF;
                        helpPtr = (uint8_t*)helpPtr + 5;
                        break;
                    case 6:
                        *(uint64_t*)helpPtr = (uint64_t)val.second[k] & 0xFFFFFFFFFFFF;
                        helpPtr = (uint8_t*)helpPtr + 6;
                        break;
                    case 7:
                        *(uint64_t*)helpPtr = (uint64_t)val.second[k] & 0xFFFFFFFFFFFFFF;
                        helpPtr = (uint8_t*)helpPtr + 7;
                        break;
                    case 8:
                        *(uint64_t*)helpPtr = (uint64_t)val.second[k];
                        helpPtr = (uint64_t*)helpPtr + 1;
                        break;
                    }
                }

                // write a delimiter of the correct width
                // switch (val.second[val.second.size() - 1]) {
                // case 1:
                //     *(uint8_t*)helpPtr = (uint8_t)DELIM;
                //     helpPtr = (uint8_t*)helpPtr + 1;
                //     break;
                // case 2:
                //     *(uint16_t*)helpPtr = (uint16_t)DELIM;
                //     helpPtr = (uint16_t*)helpPtr + 1;
                //     break;
                // case 4:
                //     *(uint32_t*)helpPtr = (uint32_t)DELIM;
                //     helpPtr = (uint32_t*)helpPtr + 1;
                //     break;
                // case 8:
                //     *(uint64_t*)helpPtr = (uint64_t)DELIM;
                //     helpPtr = (uint64_t*)helpPtr + 1;
                //     break;
                // }
                memset(helpPtr, 0, val.second[val.second.size() - 1]);
                helpPtr = (uint8_t*)helpPtr + val.second[val.second.size() - 1];
            }
        }

        //* end logic

        metadata = new uint32_t[NUM_META_DATA];

        // Set the meta data
        metadata[0] = compressionLevel;
        metadata[1] = innerDim;
        metadata[2] = outerDim;
        metadata[3] = nnz;
        metadata[4] = val_t;
        metadata[5] = index_t;

        // calculate the compression size
        calculateCompSize();

        // run the user checks
        #ifdef IVSPARSE_DEBUG
        userChecks();
        #endif

    }  // end of private transpose constructor

}  // namespace IVSparse