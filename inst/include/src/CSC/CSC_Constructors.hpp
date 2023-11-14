/**
 * @file CSC_Constructors.hpp
 * @author Skyler Ruiter and Seth Wolfgang
 * @brief Constructors for CSC Sparse Matrices
 * @version 0.1
 * @date 2023-07-03
 */

#pragma once

namespace IVSparse {

// Destructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::~SparseMatrix() {
  if (vals != nullptr) {
    free(vals);
  }
  if (innerIdx != nullptr) {
    free(innerIdx);
  }
  if (outerPtr != nullptr) {
    free(outerPtr);
  }
  if (metadata != nullptr) {
    delete[] metadata;
  }
}

// Eigen Constructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::SparseMatrix(Eigen::SparseMatrix<T> &mat) {
  
  // Make sure the matrix is compressed before reading it in
  mat.makeCompressed();

  // set the dimensions and nnz
  innerDim = mat.innerSize();
  outerDim = mat.outerSize();
  numRows = mat.rows();
  numCols = mat.cols();
  nnz = mat.nonZeros();

  encodeValueType();
  index_t = sizeof(indexT);

  // allocate the memory
  try {
    vals = (T *)malloc(nnz * sizeof(T));
    innerIdx = (indexT *)malloc(nnz * sizeof(indexT));
    outerPtr = (indexT *)malloc((outerDim + 1) * sizeof(indexT));
  } catch (std::bad_alloc &e) {
    std::cerr << "Allocation failed: " << e.what() << '\n';
  }

  // set the metadata
  metadata = new uint32_t[NUM_META_DATA];
  metadata[0] = 1;
  metadata[1] = innerDim;
  metadata[2] = outerDim;
  metadata[3] = nnz;
  metadata[4] = val_t;
  metadata[5] = index_t;

  // copy the data
  memcpy(vals, mat.valuePtr(), sizeof(T) * nnz);
  memcpy(innerIdx, mat.innerIndexPtr(), sizeof(indexT) * nnz);
  memcpy(outerPtr, mat.outerIndexPtr(), sizeof(indexT) * (outerDim + 1));

  // calculate the compressed size and run the user checks
  calculateCompSize();

  // run the user checks
  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif
}

// eigen sparse matrix constructor (row major)
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::SparseMatrix(Eigen::SparseMatrix<T, Eigen::RowMajor> &other) {
  
  other.makeCompressed();

  // set the dimensions and nnz
  innerDim = other.innerSize();
  outerDim = other.outerSize();
  numRows = other.rows();
  numCols = other.cols();
  nnz = other.nonZeros();

  encodeValueType();
  index_t = sizeof(indexT);

  // allocate the memory
  try {
    vals = (T *)malloc(nnz * sizeof(T));
    innerIdx = (indexT *)malloc(nnz * sizeof(indexT));
    outerPtr = (indexT *)malloc((outerDim + 1) * sizeof(indexT));
  } catch (std::bad_alloc &e) {
    std::cerr << "Allocation failed: " << e.what() << '\n';
  }

  // set the metadata
  metadata = new uint32_t[NUM_META_DATA];
  metadata[0] = 1;
  metadata[1] = innerDim;
  metadata[2] = outerDim;
  metadata[3] = nnz;
  metadata[4] = val_t;
  metadata[5] = index_t;

  // copy the data
  memcpy(vals, other.valuePtr(), sizeof(T) * nnz);
  memcpy(innerIdx, other.innerIndexPtr(), sizeof(indexT) * nnz);
  memcpy(outerPtr, other.outerIndexPtr(), sizeof(indexT) * (outerDim + 1));

  // calculate the compressed size and run the user checks
  calculateCompSize();

  // run the user checks
  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif
}

// Deep Copy Constructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::SparseMatrix(const IVSparse::SparseMatrix<T, indexT, 1, columnMajor> &other) {
  *this = other;
}

// General Conversion Constructor
template <typename T, typename indexT, bool columnMajor>
template <uint8_t compressionLevel2>
SparseMatrix<T, indexT, 1, columnMajor>::SparseMatrix(IVSparse::SparseMatrix<T, indexT, compressionLevel2, columnMajor> &other) {
  
  // if already CSC then just copy
  if constexpr (compressionLevel2 == 1) {
    *this = other;
    return;
  }

  // make a temporary CSC matrix
  IVSparse::SparseMatrix<T, indexT, 1, columnMajor> temp;

  if constexpr (compressionLevel2 == 2) {
    temp = other.toCSC();
  } else if constexpr (compressionLevel2 == 3) {
    temp = other.toCSC();
  }

  // copy the temporary matrix
  *this = temp;

  // run the user checks
  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif
}

// Raw CSC Constructor
template <typename T, typename indexT, bool columnMajor>
template <typename T2, typename indexT2>
SparseMatrix<T, indexT, 1, columnMajor>::SparseMatrix(
T2 *vals, indexT2 *innerIndices, indexT2 *outerPtr, uint32_t num_rows, uint32_t num_cols, uint32_t nnz) {
  
  #ifdef IVSPARSE_DEBUG
  assert(vals != nullptr && "Values array is null");
  assert(innerIndices != nullptr && "Inner indices array is null");
  assert(outerPtr != nullptr && "Outer pointer array is null");
  assert(num_rows > 0 && "Number of rows must be greater than 0");
  assert(num_cols > 0 && "Number of columns must be greater than 0");
  assert(nnz > 0 && "Number of nonzeros must be greater than 0");
  #endif

  // set the dimensions and nnz
  if constexpr (columnMajor) {
    innerDim = num_rows;
    outerDim = num_cols;
  } else {
    innerDim = num_cols;
    outerDim = num_rows;
  }

  numRows = num_rows;
  numCols = num_cols;
  this->nnz = nnz;

  encodeValueType();
  index_t = sizeof(indexT);

  // allocate the memory
  try {
    this->vals = (T *)malloc(nnz * sizeof(T));
    innerIdx = (indexT *)malloc(nnz * sizeof(indexT));
    this->outerPtr = (indexT *)malloc((outerDim + 1) * sizeof(indexT));
  } catch (std::bad_alloc &e) {
    std::cerr << "Allocation failed: " << e.what() << '\n';
  }

  // set the metadata
  metadata = new uint32_t[NUM_META_DATA];
  metadata[0] = 1;
  metadata[1] = innerDim;
  metadata[2] = outerDim;
  metadata[3] = nnz;
  metadata[4] = val_t;
  metadata[5] = index_t;

  // copy the data
  memcpy(this->vals, vals, sizeof(T) * nnz);
  memcpy(innerIdx, innerIndices, sizeof(indexT) * nnz);
  memcpy(this->outerPtr, outerPtr, sizeof(indexT) * (outerDim + 1));

  // calculate the compressed size and run the user checks
  calculateCompSize();

  // run the user checks
  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif
}

// Vector Constructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::SparseMatrix(typename IVSparse::SparseMatrix<T, indexT, 1, columnMajor>::Vector &vec) {
  
  // set the dimensions and nnz
  if constexpr (columnMajor) {
    innerDim = vec.getLength();
    outerDim = 1;
    numRows = vec.getLength();
    numCols = 1;
  } else {
    innerDim = 1;
    outerDim = vec.getLength();
    numRows = 1;
    numCols = vec.getLength();
  }

  nnz = vec.nonZeros();

  encodeValueType();
  index_t = sizeof(indexT);

  metadata = new uint32_t[NUM_META_DATA];
  metadata[0] = 1;
  metadata[1] = innerDim;
  metadata[2] = outerDim;
  metadata[3] = nnz;
  metadata[4] = val_t;
  metadata[5] = index_t;

  // if the vector is empty, return
  if (vec.byteSize() == 0) {
    vals = nullptr;
    innerIdx = nullptr;

    // allocate outerPtr
    try {
      outerPtr = (indexT *)malloc((outerDim + 1) * sizeof(indexT));
    } catch (std::bad_alloc &e) {
      std::cerr << "Allocation failed: " << e.what() << '\n';
    }

    outerPtr[0] = 0;
    return;
  }

  // allocate the memory
  try {
    vals = (T *)malloc(nnz * sizeof(T));
    innerIdx = (indexT *)malloc(nnz * sizeof(indexT));
    outerPtr = (indexT *)malloc((outerDim + 1) * sizeof(indexT));
  } catch (std::bad_alloc &e) {
    std::cerr << "Allocation failed: " << e.what() << '\n';
  }

  // copy the data and set the outerPtr
  memcpy(vals, vec.getValues(), sizeof(T) * nnz);
  memcpy(innerIdx, vec.getInnerIndices(), sizeof(indexT) * nnz);
  outerPtr[0] = 0;
  outerPtr[1] = nnz;

  // calculate the compressed size and run the user checks
  calculateCompSize();

  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif
}

// Array of Vectors Constructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::SparseMatrix(
std::vector<typename IVSparse::SparseMatrix<T, indexT, 1, columnMajor>::Vector> &vecs) {
  
  // Make a temporary CSC matrix
  IVSparse::SparseMatrix<T, indexT, 1, columnMajor> temp(vecs[0]);

  // append on each vector in the array
  for (size_t i = 1; i < vecs.size(); i++) { temp.append(vecs[i]); }

  // copy the temporary matrix
  *this = temp;

  // run the user checks and calculate the compressed size
  calculateCompSize();

  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif
}

// File Constructor
template <typename T, typename indexT, bool columnMajor>
SparseMatrix<T, indexT, 1, columnMajor>::SparseMatrix(const char *filename) {
  
  FILE *fp = fopen(filename, "rb");

  #ifdef IVSPARSE_DEBUG
  // make sure the file exists
  if (fp == NULL) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    exit(1);
  }
  #endif

  // read the metadata and set the dimensions/nnz
  metadata = new uint32_t[NUM_META_DATA];
  fread(metadata, sizeof(uint32_t), NUM_META_DATA, fp);
  innerDim = metadata[1];
  outerDim = metadata[2];
  nnz = metadata[3];
  val_t = metadata[4];
  index_t = metadata[5];

  if constexpr (columnMajor) {
    numRows = innerDim;
    numCols = outerDim;
  } else {
    numRows = outerDim;
    numCols = innerDim;
  }

  // allocate the memory
  try {
    vals = (T *)malloc(nnz * sizeof(T));
    innerIdx = (indexT *)malloc(nnz * sizeof(indexT));
    outerPtr = (indexT *)malloc((outerDim + 1) * sizeof(indexT));
  } catch (std::bad_alloc &e) {
    std::cerr << "Error: Could not allocate memory for IVSparse matrix"
              << std::endl;
    exit(1);
  }

  // read the data
  fread(vals, sizeof(T), nnz, fp);
  fread(innerIdx, sizeof(indexT), nnz, fp);
  fread(outerPtr, sizeof(indexT), outerDim + 1, fp);

  // close the file
  fclose(fp);

  // run the user checks
  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif  

  calculateCompSize();
}

// COO -> CSC constructor
template <typename T, typename indexT, bool columnMajor>
template <typename T2, typename indexT2>
SparseMatrix<T, indexT, 1, columnMajor>::SparseMatrix(
std::vector<std::tuple<indexT2, indexT2, T2>> &entries, uint32_t num_rows, uint32_t num_cols, uint32_t nnz) {
  
  #ifdef IVSPARSE_DEBUG
  assert(entries.size() > 0 && "Entries array is empty");
  assert(num_rows > 0 && "Number of rows must be greater than 0");
  assert(num_cols > 0 && "Number of columns must be greater than 0");
  assert(nnz > 0 && "Number of nonzeros must be greater than 0");
  #endif

  if constexpr (columnMajor) {
    innerDim = num_rows;
    outerDim = num_cols;
  } else {
    innerDim = num_cols;
    outerDim = num_rows;
  }

  numRows = num_rows;
  numCols = num_cols;

  this->nnz = nnz;

  encodeValueType();
  index_t = sizeof(indexT);

  try {
    vals = (T *)malloc(nnz * sizeof(T));
    innerIdx = (indexT *)malloc(nnz * sizeof(indexT));
    outerPtr = (indexT *)calloc(outerDim + 1, sizeof(indexT));
  } catch (std::bad_alloc &e) {
    std::cerr << "Allocation failed: " << e.what() << '\n';
  }

  metadata = new uint32_t[NUM_META_DATA];

  metadata[0] = 1;
  metadata[1] = innerDim;
  metadata[2] = outerDim;
  metadata[3] = nnz;
  metadata[4] = val_t;
  metadata[5] = index_t;

  // sort the tuples by first by column then by row
  std::sort(entries.begin(), entries.end(),
            [](const std::tuple<indexT2, indexT2, T2> &a,
               const std::tuple<indexT2, indexT2, T2> &b) {
              if (std::get<1>(a) == std::get<1>(b)) {
                return std::get<0>(a) < std::get<0>(b);
              } else {
                return std::get<1>(a) < std::get<1>(b);
              }
            });

  int OuterIndex = -1;
  int count = 0;

  for (size_t i = 0; i < entries.size(); i++) {
    vals[i] = std::get<2>(entries[i]);
    innerIdx[i] = std::get<0>(entries[i]);
    count++;

    if (std::get<1>(entries[i]) != OuterIndex) {
      if (OuterIndex != -1) {
        outerPtr[OuterIndex + 1] = count - 1;
      }
      OuterIndex = std::get<1>(entries[i]);
      outerPtr[OuterIndex] = count - 1;
    }
  }

  outerPtr[OuterIndex + 1] = count;

  // run the user checks
  #ifdef IVSPARSE_DEBUG
  userChecks();
  #endif

  calculateCompSize();
}

}  // namespace IVSparse