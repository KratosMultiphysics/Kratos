---
title: CsrMatrix
keywords: CsrMatrix, Sparse Matrix, CSR
tags: [CsrMatrix, Data Structures, Linear Algebra]
sidebar: kratos_for_developers
summary: Documentation for the CsrMatrix class, a serial Compressed Sparse Row matrix implementation in Kratos.
---

## CsrMatrix Class

The `CsrMatrix` class in Kratos implements a serial (non-distributed) Compressed Sparse Row (CSR) matrix. CSR is an efficient format for storing sparse matrices, reducing memory consumption and allowing for fast matrix-vector operations. This class is templated, allowing for different data types for its values and indices.

### Template Parameters

*   `TDataType`: The data type of the values stored in the matrix (e.g., `double`, `float`). Defaults to `double`.
*   `TIndexType`: The data type used for row and column indices (e.g., `std::size_t`, `int`). Defaults to `std::size_t`.

### Lifecycle Methods

The `CsrMatrix` can be created in several ways, depending on whether the sparsity pattern is known beforehand and how the data is structured.

#### Constructors

*   `CsrMatrix()`: Default constructor. Creates an empty matrix associated with the serial `DataCommunicator`.
    ```cpp
    // KratosCoreFastSuite :: CSRMatrixConstruction
    // (Simplified for clarity)
    Kratos::CsrMatrix<double> empty_matrix; 
    KRATOS_EXPECT_EQ(empty_matrix.size1(), 0); // Example check
    KRATOS_EXPECT_EQ(empty_matrix.size2(), 0); // Example check
    KRATOS_EXPECT_EQ(empty_matrix.nnz(), 0);   // Example check
    ```

*   `CsrMatrix(const DataCommunicator& rComm)`: Constructor that takes a `DataCommunicator`. It will raise an error if a distributed communicator is passed, as `CsrMatrix` is serial.
    ```cpp
    // Assuming Kratos::ParallelEnvironment::GetDefaultDataCommunicator() is serial
    const Kratos::DataCommunicator& serial_comm = Kratos::ParallelEnvironment::GetDefaultDataCommunicator();
    Kratos::CsrMatrix<double> matrix_with_comm(serial_comm);
    ```

*   `template<class TGraphType> CsrMatrix(const TGraphType& rSparseGraph)`: Constructs a `CsrMatrix` from a given sparse graph (e.g., `SparseContiguousRowGraph` or `SparseGraph`). The matrix structure (row and column indices) is copied from the graph, and values are initialized to zero. The communicator is also taken from the graph. This is a common way to initialize a matrix when its sparsity pattern is known from a graph.
    ```cpp
    // KratosCoreFastSuite :: CSRMatrixConstruction (Conceptual Example)
    // Assume 'connectivities' is a std::vector of std::vector<IndexType>
    // representing element connectivities.
    Kratos::SparseContiguousRowGraph<Kratos::Globals::DataLocation::Local> Agraph(40); // Number of nodes (example value)

    // Populate Agraph based on connectivities (e.g., from a mesh)
    // This part is highly dependent on how 'connectivities' are defined and added.
    // For example, if connectivities = {{0,1}, {1,2}, {0,2}} for a 3-node graph:
    std::vector<std::vector<Kratos::CsrMatrix<double>::IndexType>> connectivities_example = {{0,1},{1,2},{0,2}};
    for(const auto& c : connectivities_example) {
        Agraph.AddEntries(c); // Add symmetric entries for each connectivity
    }
    Agraph.Finalize(); // Finalize graph structure

    // Create CsrMatrix from the graph
    Kratos::CsrMatrix<double> A(Agraph); 
    // Matrix A now has the sparsity pattern of Agraph, with all values set to 0.0.
    // A.size1() would be 40 (or the number of nodes used to initialize Agraph).
    // A.nnz() would be the number of unique non-zero locations defined by Agraph.
    ```
    *Note: In the actual test `CSRMatrixConstruction`, `Agraph` is populated from connectivities and then used to construct `A`. Values are then assembled into `A`.*

*   `CsrMatrix(const MatrixMapType& rMatrixMap)`: Constructs a `CsrMatrix` from a `MatrixMapType`. `MatrixMapType` is an alias for `std::unordered_map<std::pair<IndexType, IndexType>, TDataType, ...>`, which maps `(row, col)` pairs to values. The structure and values are determined from the map.
    ```cpp
    // Conceptual Example
    Kratos::CsrMatrix<double>::MatrixMapType matrix_data_map;
    matrix_data_map[{0, 0}] = 1.0;
    matrix_data_map[{0, 1}] = 2.0;
    matrix_data_map[{1, 1}] = 3.0;

    Kratos::CsrMatrix<double> matrix_from_map(matrix_data_map);
    // matrix_from_map will have:
    // (0,0) = 1.0
    // (0,1) = 2.0
    // (1,1) = 3.0
    // size1() and size2() will be determined by the max indices in the map.
    // For this example, size1() would likely be 2, size2() would be 2.
    KRATOS_EXPECT_EQ(matrix_from_map(0,0), 1.0); // Example check
    KRATOS_EXPECT_EQ(matrix_from_map(0,1), 2.0); // Example check
    KRATOS_EXPECT_EQ(matrix_from_map(1,1), 3.0); // Example check
    ```

*   `CsrMatrix(const CsrMatrix<TDataType,TIndexType>& rOtherMatrix)`: Copy constructor. Performs a deep copy of the other matrix's data.
    ```cpp
    Kratos::CsrMatrix<double> original_matrix; 
    // ... (original_matrix is initialized and filled) ...
    // For example:
    Kratos::CsrMatrix<double>::MatrixMapType temp_map;
    temp_map[{0,0}] = 5.0; temp_map[{1,1}] = 6.0;
    original_matrix = Kratos::CsrMatrix<double>(temp_map);

    Kratos::CsrMatrix<double> copied_matrix(original_matrix);
    // copied_matrix now has its own copy of the data from original_matrix.
    KRATOS_EXPECT_EQ(copied_matrix(0,0), original_matrix(0,0)); // Example check
    ```

*   `CsrMatrix(CsrMatrix<TDataType,TIndexType>&& rOtherMatrix)`: Move constructor. Takes ownership of the data from `rOtherMatrix`. `rOtherMatrix` is left in a valid but unspecified state (typically empty).
    ```cpp
    Kratos::CsrMatrix<double> source_matrix;
    // ... (source_matrix is initialized and filled) ...
    Kratos::CsrMatrix<double>::MatrixMapType temp_map_move;
    temp_map_move[{0,0}] = 7.0;
    source_matrix = Kratos::CsrMatrix<double>(temp_map_move);
    Kratos::CsrMatrix<double>::IndexType original_nnz = source_matrix.nnz();

    Kratos::CsrMatrix<double> moved_matrix(std::move(source_matrix));
    // moved_matrix now owns the data previously in source_matrix.
    KRATOS_EXPECT_EQ(moved_matrix.nnz(), original_nnz); // Example check
    KRATOS_EXPECT_EQ(source_matrix.nnz(), 0);           // Example check: source is empty
    // source_matrix should not be used unless reinitialized.
    ```

#### Destructor

*   `~CsrMatrix()`: Destructor. Frees the allocated memory for row indices, column indices, and values if the matrix owns its data (i.e., `IsOwnerOfData()` is `true`).

#### Assignment Operators

*   `CsrMatrix& operator=(CsrMatrix const& rOtherMatrix) = delete;`: The copy assignment operator is deleted to prevent accidental expensive copying. If a copy is needed, use the copy constructor explicitly.

*   `CsrMatrix& operator=(CsrMatrix&& rOtherMatrix)`: Move assignment operator. Takes ownership of the data from `rOtherMatrix`.
    ```cpp
    Kratos::CsrMatrix<double> matrix_to_assign_to;
    Kratos::CsrMatrix<double> another_matrix;
    // ... (another_matrix is initialized and filled) ...
    Kratos::CsrMatrix<double>::MatrixMapType temp_map_assign;
    temp_map_assign[{0,1}] = 8.0;
    another_matrix = Kratos::CsrMatrix<double>(temp_map_assign);
    Kratos::CsrMatrix<double>::IndexType original_nnz_assign = another_matrix.nnz();

    matrix_to_assign_to = std::move(another_matrix);
    // matrix_to_assign_to now owns the data from another_matrix.
    KRATOS_EXPECT_EQ(matrix_to_assign_to.nnz(), original_nnz_assign); // Example check
    KRATOS_EXPECT_EQ(another_matrix.nnz(), 0);                        // Example check
    ```

### Main Operators and Member Functions

This section details various functions for interacting with and manipulating the `CsrMatrix`.

#### General Utility

*   `void Clear()`: Clears all data in the matrix, resetting it to an empty state (0 rows, 0 columns, no non-zero entries). Deallocates memory if owned.
    ```cpp
    Kratos::CsrMatrix<double> matrix_to_clear;
    // ... (matrix_to_clear is initialized and filled) ...
    Kratos::CsrMatrix<double>::MatrixMapType temp_map_clear;
    temp_map_clear[{0,0}] = 1.0;
    matrix_to_clear = Kratos::CsrMatrix<double>(temp_map_clear);

    matrix_to_clear.Clear();
    KRATOS_EXPECT_EQ(matrix_to_clear.size1(), 0);
    KRATOS_EXPECT_EQ(matrix_to_clear.nnz(), 0);
    ```

*   `const DataCommunicator& GetComm() const`: Returns a constant reference to the `DataCommunicator` associated with the matrix.
*   `const DataCommunicator* pGetComm() const`: Returns a pointer to the `DataCommunicator` associated with the matrix.

*   `void SetValue(const TDataType value)`: Sets all *existing* non-zero entries (i.e., entries present in the sparsity pattern) in the matrix to the specified `value`. This does not change the sparsity pattern.
    ```cpp
    // Assume 'A' is a CsrMatrix initialized with a sparsity pattern
    Kratos::CsrMatrix<double>::MatrixMapType sv_map;
    sv_map[{0,0}]=1.0; sv_map[{1,1}]=1.0;
    Kratos::CsrMatrix<double> A_setval(sv_map);
    A_setval.SetValue(10.0); // All stored entries in A_setval are now 10.0
    KRATOS_EXPECT_EQ(A_setval(0,0), 10.0); // Example check
    KRATOS_EXPECT_EQ(A_setval(1,1), 10.0); // Example check
    ```

*   `IndexType size1() const`: Returns the number of rows in the matrix.
    ```cpp
    // Given CsrMatrix A
    Kratos::CsrMatrix<double> A_s1; // Assume A_s1 is 3xN
    Kratos::CsrMatrix<double>::MatrixMapType s1_map;
    s1_map[{0,0}]=1; s1_map[{1,0}]=1; s1_map[{2,0}]=1;
    A_s1 = Kratos::CsrMatrix<double>(s1_map);
    Kratos::CsrMatrix<double>::IndexType num_rows = A_s1.size1(); // num_rows would be 3
    ```

*   `IndexType size2() const`: Returns the number of columns in the matrix.
    ```cpp
    // Given CsrMatrix A
    Kratos::CsrMatrix<double> A_s2; // Assume A_s2 is Mx2
    Kratos::CsrMatrix<double>::MatrixMapType s2_map;
    s2_map[{0,0}]=1; s2_map[{0,1}]=1;
    A_s2 = Kratos::CsrMatrix<double>(s2_map);
    Kratos::CsrMatrix<double>::IndexType num_cols = A_s2.size2(); // num_cols would be 2
    ```
    *Note: `size2()` is typically set via `SetColSize()` or `ComputeColSize()` after the sparsity pattern is defined.*

*   `IndexType nnz() const`: Returns the number of non-zero entries (or stored entries, as CSR can store explicit zeros if they are part of the pattern).
    ```cpp
    // Given CsrMatrix A
    Kratos::CsrMatrix<double>::MatrixMapType nnz_map;
    nnz_map[{0,0}]=1; nnz_map[{0,1}]=2; nnz_map[{1,0}]=3;
    Kratos::CsrMatrix<double> A_nnz(nnz_map);
    Kratos::CsrMatrix<double>::IndexType non_zeros = A_nnz.nnz(); // non_zeros would be 3
    ```

*   `bool IsOwnerOfData() const`: Returns `true` if the matrix owns the underlying data arrays (row indices, column indices, values), `false` otherwise.
*   `void SetIsOwnerOfData(bool IsOwner)`: Sets whether the matrix owns its data. Use with extreme caution, as this can lead to memory issues if not managed correctly. Primarily for advanced use cases where matrix data might be shared or managed externally.

#### Data Access (CSR Arrays)

These methods provide access to the raw CSR data arrays. Modifying these directly requires a deep understanding of the CSR format.

*   `Kratos::span<IndexType>& index1_data()`: Returns a writable `span` for the row pointers (CSR `row_ptr` array). `index1_data()[i]` is the start index in `index2_data()` and `value_data()` for row `i`. The size of this span is `size1() + 1`.
*   `const Kratos::span<IndexType>& index1_data() const`: Const version.
*   `Kratos::span<IndexType>& index2_data()`: Returns a writable `span` for the column indices (CSR `col_ind` array).
*   `const Kratos::span<IndexType>& index2_data() const`: Const version.
*   `Kratos::span<TDataType>& value_data()`: Returns a writable `span` for the values of the stored entries.
*   `const Kratos::span<TDataType>& value_data() const`: Const version.

    ```cpp
    // Example: Iterating through entries of row `i` (conceptual)
    // Assume A_iter is an initialized CsrMatrix
    Kratos::CsrMatrix<double>::MatrixMapType iter_map;
    iter_map[{0,0}]=1.0; iter_map[{0,1}]=2.0; iter_map[{1,1}]=3.0;
    Kratos::CsrMatrix<double> A_iter(iter_map);
    Kratos::CsrMatrix<double>::IndexType row_i = 0; // some row index
    if (row_i < A_iter.size1()) {
        const auto& row_indices = A_iter.index1_data();
        const auto& col_indices = A_iter.index2_data();
        const auto& values = A_iter.value_data();

        Kratos::CsrMatrix<double>::IndexType row_start = row_indices[row_i];
        Kratos::CsrMatrix<double>::IndexType row_end = row_indices[row_i + 1];

        for (Kratos::CsrMatrix<double>::IndexType k = row_start; k < row_end; ++k) {
            Kratos::CsrMatrix<double>::IndexType col_j = col_indices[k];
            Kratos::CsrMatrix<double>::DataType value = values[k];
            // Process (row_i, col_j, value)
            std::cout << "A(" << row_i << "," << col_j << ") = " << value << std::endl; // Example action
        }
    }
    ```

#### Sizing and Capacity

*   `void SetColSize(IndexType Ncols)`: Sets the number of columns (`mNcols`) reported by `size2()`. Does not allocate or change data arrays.
*   `void SetRowSize(IndexType Nrows)`: Sets the number of rows (`mNrows`) reported by `size1()`. Does not allocate or change data arrays, but `index1_data` should be consistent (size `Nrows+1`).
*   `void ComputeColSize()`: Calculates `mNcols` (the value returned by `size2()`) by finding the maximum column index present in `index2_data()` and adding 1. Call this if you manually populate `index2_data()` or construct the matrix from external CSR arrays without initially setting `mNcols`.
    ```cpp
    // Assume matrix A_compcol is constructed, and its column indices are populated.
    Kratos::CsrMatrix<double>::MatrixMapType compcol_map;
    compcol_map[{0,0}]=1; compcol_map[{0,2}]=1; // Max column index is 2
    Kratos::CsrMatrix<double> A_compcol(compcol_map);
    A_compcol.ComputeColSize(); 
    // Now A_compcol.size2() will reflect the maximum column index found + 1 (i.e., 3).
    KRATOS_EXPECT_EQ(A_compcol.size2(), 3); // Example check
    ```
*   `void CheckColSize()`: Verifies if any column index in `mColIndices` is greater than or equal to `mNcols`. Throws an error if an inconsistency is found. Useful for debugging.

#### Memory Management (Raw Data Pointers)

These methods are for advanced use cases, allowing the `CsrMatrix` to use externally managed memory or to resize its own memory buffers. **Use with extreme caution.** If `IsOwnerOfData()` is `false`, the user is responsible for managing the lifetime of these external arrays. If `IsOwnerOfData()` is `true`, resizing will deallocate old data and allocate new data.

*   `void AssignIndex1Data(TIndexType* pExternalData, TIndexType DataSize)`: Assigns an external array `pExternalData` (e.g., from a C-style array or other memory source) of size `DataSize` as the row pointers array. The matrix will not own this data.
*   `void AssignIndex2Data(TIndexType* pExternalData, TIndexType DataSize)`: Assigns an external array for column indices.
*   `void AssignValueData(TDataType* pExternalData, TIndexType DataSize)`: Assigns an external array for values.
    ```cpp
    // Conceptual: Using externally allocated data
    Kratos::CsrMatrix<double>::IndexType num_rows_ext = 2;
    Kratos::CsrMatrix<double>::IndexType nnz_count_ext = 3;

    Kratos::CsrMatrix<double>::IndexType* my_row_ptr = new Kratos::CsrMatrix<double>::IndexType[num_rows_ext + 1]{0, 2, 3};
    Kratos::CsrMatrix<double>::IndexType* my_col_ind = new Kratos::CsrMatrix<double>::IndexType[nnz_count_ext]{0, 1, 1};
    double* my_values = new double[nnz_count_ext]{1.0, 2.0, 3.0};

    Kratos::CsrMatrix<double> A_ext;
    A_ext.SetIsOwnerOfData(false); // Crucial!
    A_ext.AssignIndex1Data(my_row_ptr, num_rows_ext + 1);
    A_ext.AssignIndex2Data(my_col_ind, nnz_count_ext);
    A_ext.AssignValueData(my_values, nnz_count_ext);
    A_ext.SetRowSize(num_rows_ext); // Manually set sizes
    A_ext.SetColSize(2);      // Manually set sizes (or ComputeColSize if indices are complex)

    // ... use A_ext ...
    KRATOS_EXPECT_EQ(A_ext(1,1), 3.0); // Example check

    // User must delete the arrays when done, as A_ext does not own them.
    delete[] my_row_ptr;
    delete[] my_col_ind;
    delete[] my_values;
    ```

*   `void ResizeIndex1Data(TIndexType DataSize)`: Resizes the internal row pointers array. Only if `IsOwnerOfData()` is `true`. Existing data is lost.
*   `void ResizeIndex2Data(TIndexType DataSize)`: Resizes internal column indices array. Only if `IsOwnerOfData()` is `true`.
*   `void ResizeValueData(TIndexType DataSize)`: Resizes internal values array. Only if `IsOwnerOfData()` is `true`.
    ```cpp
    // Conceptual: Resizing owned data (e.g., before manual population)
    Kratos::CsrMatrix<double> B_resize; // Owns its data by default
    B_resize.ResizeIndex1Data(101); // For 100 rows
    B_resize.ResizeIndex2Data(500); // For up to 500 non-zeros
    B_resize.ResizeValueData(500);
    B_resize.SetRowSize(100);
    // Now B_resize.index1_data(), B_resize.index2_data(), B_resize.value_data() can be populated.
    // Remember to set B_resize.index1_data()[0] = 0; and fill row pointers correctly.
    // Then call B_resize.ComputeColSize();
    KRATOS_EXPECT_EQ(B_resize.index1_data().size(), 101); // Example check
    ```

#### Element Access

*   `TDataType& operator()(IndexType I, IndexType J)`: Returns a reference to the element at row `I` and column `J`. If the element is not found in the sparsity pattern, it **throws an error in debug mode** (or may lead to undefined behavior in release). This operation involves a binary search for column `J` within row `I`'s entries.
*   `const TDataType& operator()(IndexType I, IndexType J) const`: Const version.
*   `bool Has(IndexType I, IndexType J) const`: Checks if an entry exists at row `I` and column `J` in the sparsity pattern. Returns `true` if it exists, `false` otherwise.

    ```cpp
    // Assume CsrMatrix A_access is initialized and assembled.
    // Example from CSRMatrixSmallRectangularMatrixMultiply test
    Kratos::CsrMatrix<double>::MatrixMapType access_map;
    access_map[{0,0}]=1.0; access_map[{0,3}]=7.0; access_map[{0,4}]=2.0;
    access_map[{1,1}]=3.0;
    access_map[{2,3}]=7.0; access_map[{2,4}]=7.0;
    Kratos::CsrMatrix<double> A_access(access_map);
    // A_access might represent:
    // [[1,0,0,7,2],
    //  [0,3,0,0,0],
    //  [0,0,0,7,7]] (Note: map doesn't predefine all zeros explicitly like a graph does)

    if (A_access.Has(0, 3)) {
        double val = A_access(0, 3); // val will be 7.0
        A_access(0, 3) = 10.0;       // Modify existing entry
        KRATOS_EXPECT_EQ(A_access(0,3), 10.0); // Example check
    }

    bool has_entry = A_access.Has(1, 0); // has_entry will be false
    // double non_existent = A_access(1,0); // This would likely cause an error in debug
    ```

#### Matrix-Vector Products

*   `template<class TInputVectorType, class TOutputVectorType> void SpMV(const TInputVectorType& rX, TOutputVectorType& rY) const`: Performs **y += A\*x**.
    *   `rX`: Input vector (e.g., `SystemVector` or any type supporting `operator()` access and `size()`).
    *   `rY`: Output vector, to which `A*rX` is added. Must be correctly sized.
*   `template<class TInputVectorType, class TOutputVectorType> void SpMV(const TDataType Alpha, const TInputVectorType& rX, const TDataType Beta, TOutputVectorType& rY) const`: Performs **y = α\*A\*x + β\*y**.

    ```cpp
    // From KratosCoreFastSuite :: CSRMatrixSpMV (simplified)
    // Assume CsrMatrix 'A_spmv' is initialized (e.g., 3x3 for this example)
    Kratos::CsrMatrix<double>::MatrixMapType spmv_map;
    spmv_map[{0,0}]=2; spmv_map[{0,1}]=-1;
    spmv_map[{1,0}]=-1; spmv_map[{1,1}]=2; spmv_map[{1,2}]=-1;
    spmv_map[{2,1}]=-1; spmv_map[{2,2}]=2;
    Kratos::CsrMatrix<double> A_spmv(spmv_map); // A common stiffness matrix type
    A_spmv.ComputeColSize(); // Ensure size2 is correct

    Kratos::Vector x_vec(A_spmv.size2()); 
    Kratos::Vector y_vec(A_spmv.size1());

    // Initialize x_vec and y_vec
    for(unsigned int i=0; i<x_vec.size(); ++i) x_vec[i] = static_cast<double>(i+1); // x = [1, 2, 3]
    for(unsigned int i=0; i<y_vec.size(); ++i) y_vec[i] = 0.0; // y = [0, 0, 0]

    A_spmv.SpMV(x_vec, y_vec); // y_vec now contains A_spmv * x_vec
    // Expected y_vec = [2*1-1*2, -1*1+2*2-1*3, -1*2+2*3] = [0, 0, 4]
    KRATOS_EXPECT_VECTOR_NEAR(y_vec, Kratos::Vector{0.0, 0.0, 4.0}, 1e-9); // Example check

    // For y = alpha*A*x + beta*y
    double alpha = 2.0;
    double beta = 0.5;
    // Re-initialize y_vec for clarity of the beta term
    for(unsigned int i=0; i<y_vec.size(); ++i) y_vec[i] = static_cast<double>(i); // y_initial = [0,1,2]
    // y_vec will be 2.0 * A_spmv * x_vec + 0.5 * y_initial
    A_spmv.SpMV(alpha, x_vec, beta, y_vec); 
    // Expected y_vec = 2.0*[0,0,4] + 0.5*[0,1,2] = [0, 0, 8] + [0, 0.5, 1] = [0, 0.5, 9]
    KRATOS_EXPECT_VECTOR_NEAR(y_vec, Kratos::Vector{0.0, 0.5, 9.0}, 1e-9); // Example check
    ```

*   `template<class TInputVectorType, class TOutputVectorType> void TransposeSpMV(const TInputVectorType& rX, TOutputVectorType& rY) const`: Performs **y += A<sup>T</sup>\*x**.
*   `template<class TInputVectorType, class TOutputVectorType> void TransposeSpMV(const TDataType Alpha, const TInputVectorType& rX, const TDataType Beta, TOutputVectorType& rY) const`: Performs **y = α\*A<sup>T</sup>\*x + β\*y**.
    *   Note: `rY` must have size `A.size2()` and `rX` must have size `A.size1()`.

    ```cpp
    // From KratosCoreFastSuite :: CSRMatrixRectangularMatrix (conceptual)
    // Assume CsrMatrix 'A_tspmv' is initialized (e.g. 2x3)
    Kratos::CsrMatrix<double>::MatrixMapType tspmv_map;
    tspmv_map[{0,0}]=1; tspmv_map[{0,1}]=2; tspmv_map[{0,2}]=3;
    tspmv_map[{1,0}]=4; tspmv_map[{1,1}]=5; tspmv_map[{1,2}]=6;
    Kratos::CsrMatrix<double> A_tspmv(tspmv_map);
    A_tspmv.ComputeColSize(); // Ensure size2 is correct (3)

    Kratos::Vector x_for_transpose(A_tspmv.size1()); // Input vector for A^T * x, size 2
    Kratos::Vector y_from_transpose(A_tspmv.size2());   // Output vector, size 3

    for(unsigned int i=0; i<x_for_transpose.size(); ++i) x_for_transpose[i] = 1.0; // x_for_transpose = [1,1]
    for(unsigned int i=0; i<y_from_transpose.size(); ++i) y_from_transpose[i] = 0.0; // y_from_transpose = [0,0,0]

    A_tspmv.TransposeSpMV(x_for_transpose, y_from_transpose); // y_from_transpose = A_tspmv^T * x_for_transpose
    // A_tspmv^T = [[1,4],[2,5],[3,6]]
    // y_from_transpose = [[1,4],[2,5],[3,6]] * [1,1]^T = [1*1+4*1, 2*1+5*1, 3*1+6*1]^T = [5,7,9]^T
    KRATOS_EXPECT_VECTOR_NEAR(y_from_transpose, Kratos::Vector{5.0,7.0,9.0}, 1e-9); // Example check
    ```

#### Norms

*   `TDataType NormFrobenius() const`: Calculates and returns the Frobenius norm of the matrix (square root of the sum of squares of its elements).
    ```cpp
    // Given CsrMatrix A_frob
    Kratos::CsrMatrix<double>::MatrixMapType frob_map;
    frob_map[{0,0}]=1; frob_map[{0,1}]=2;
    frob_map[{1,0}]=3; frob_map[{1,1}]=4;
    Kratos::CsrMatrix<double> A_frob(frob_map);
    double frob_norm = A_frob.NormFrobenius();
    // Expected: sqrt(1^2 + 2^2 + 3^2 + 4^2) = sqrt(1+4+9+16) = sqrt(30) approx 5.477
    KRATOS_EXPECT_NEAR(frob_norm, std::sqrt(30.0), 1e-9); // Example check
    ```

### Assembly Operations

These methods are crucial for constructing the matrix, especially in Finite Element Method (FEM) applications where a global matrix is built by summing contributions from local element matrices. The matrix's sparsity pattern (i.e., which entries `(i,j)` can be non-zero) MUST be defined before assembly. This is typically done by constructing the `CsrMatrix` from a `SparseGraph` or `SparseContiguousRowGraph`.

*   `void reserve(IndexType NRows, IndexType nnz_prediction)`: Pre-allocates memory if you intend to build the CSR structure more manually (e.g., by directly populating `index1_data`, `index2_data`, `value_data`).
    *   `NRows`: The number of rows the matrix will have. `index1_data` will be sized to `NRows + 1`.
    *   `nnz_prediction`: An estimate of the number of non-zero entries. `index2_data` and `value_data` will be sized to this.
    *   If `IsOwnerOfData()` is true (default), this allocates new arrays.
    ```cpp
    Kratos::CsrMatrix<double> A_manual_reserve;
    Kratos::CsrMatrix<double>::IndexType num_rows_reserve = 3;
    Kratos::CsrMatrix<double>::IndexType estimated_nnz_reserve = 5;
    A_manual_reserve.reserve(num_rows_reserve, estimated_nnz_reserve);
    A_manual_reserve.SetRowSize(num_rows_reserve);
    // A_manual_reserve.index1_data() is size 4 (num_rows_reserve + 1)
    // A_manual_reserve.index2_data() is size 5
    // A_manual_reserve.value_data() is size 5
    // Now, one would manually populate these arrays and then call ComputeColSize().
    KRATOS_EXPECT_EQ(A_manual_reserve.index1_data().size(), num_rows_reserve + 1); // Example check
    ```

*   `MatrixMapType ToMap() const`: Converts the `CsrMatrix` to a `MatrixMapType` (an `std::unordered_map<std::pair<IndexType, IndexType>, TDataType, ...>`). This can be useful for debugging, converting to other formats, or if you need easy lookup of arbitrary elements (though it's less efficient than direct CSR access for iteration).
    ```cpp
    // Assume CsrMatrix 'A_tomap' is initialized and assembled.
    Kratos::CsrMatrix<double>::MatrixMapType tomap_orig_map;
    tomap_orig_map[{0,0}]=1.0; tomap_orig_map[{1,1}]=2.0;
    Kratos::CsrMatrix<double> A_tomap(tomap_orig_map);

    Kratos::CsrMatrix<double>::MatrixMapType A_map_result = A_tomap.ToMap();
    for(const auto& entry : A_map_result) {
        Kratos::CsrMatrix<double>::IndexType r = entry.first.first;
        Kratos::CsrMatrix<double>::IndexType c = entry.first.second;
        double val = entry.second;
        std::cout << "(" << r << ", " << c << "): " << val << std::endl;
    }
    KRATOS_EXPECT_EQ(A_map_result.size(), 2); // Example check
    ```

*   `void BeginAssemble()`: In the serial `CsrMatrix` implementation, this function currently does nothing. It's a placeholder for compatibility with potential MPI (distributed memory) versions where it might initiate communication or locking mechanisms required before assembly.
*   `void FinalizeAssemble()`: Similar to `BeginAssemble()`, this function does nothing in the serial version. It's for MPI compatibility.

    ```cpp
    // Typical assembly workflow:
    Kratos::SparseContiguousRowGraph<> Agraph_assembly(3); // Example graph
    Agraph_assembly.AddEntry(0,0); Agraph_assembly.AddEntry(1,1); Agraph_assembly.AddEntry(2,2);
    Agraph_assembly.Finalize();
    Kratos::CsrMatrix<double> A_workflow(Agraph_assembly); // Assume Agraph_assembly defines the sparsity
    A_workflow.BeginAssemble(); // No-op in serial, but good practice
    // ... calls to Assemble() or AssembleEntry() ...
    A_workflow.AssembleEntry(1.0,0,0);
    A_workflow.FinalizeAssemble(); // No-op in serial
    ```

*   `template<class TMatrixType, class TIndexVectorType> void Assemble(const TMatrixType& rMatrixInput, const TIndexVectorType& rEquationId)`: Assembles (adds) a local (element) matrix `rMatrixInput` into the global `CsrMatrix`.
    *   `rMatrixInput`: The local matrix (e.g., a `Kratos::Matrix` from an element calculation). Assumed to be square.
    *   `rEquationId`: A vector mapping local indices of `rMatrixInput` (from `0` to `rMatrixInput.size1()-1`) to the global row AND column indices in `*this`.
    *   Values from `rMatrixInput` are atomically added to the corresponding entries in `*this`. The entries must exist in the pre-defined sparsity pattern.

    ```cpp
    // From KratosCoreFastSuite :: CSRMatrixConstruction (conceptual)
    // Assume CsrMatrix 'A_glob_assem' is constructed from 'Agraph_glob_assem'.
    Kratos::SparseContiguousRowGraph<> Agraph_glob_assem(3);
    std::vector<Kratos::CsrMatrix<double>::IndexType> conn_example_1 = {0,1};
    std::vector<Kratos::CsrMatrix<double>::IndexType> conn_example_2 = {1,2};
    Agraph_glob_assem.AddEntries(conn_example_1);
    Agraph_glob_assem.AddEntries(conn_example_2);
    Agraph_glob_assem.Finalize();
    Kratos::CsrMatrix<double> A_glob_assem(Agraph_glob_assem);

    // Assume 'connectivities_glob_assem' is std::vector<std::vector<IndexType>> where each
    // inner vector holds the global node IDs for an element.
    std::vector<std::vector<Kratos::CsrMatrix<double>::IndexType>> connectivities_glob_assem = {conn_example_1, conn_example_2};

    A_glob_assem.BeginAssemble();
    for(const auto& c : connectivities_glob_assem) { // 'c' is like an rEquationId vector
        Kratos::Matrix local_element_matrix(c.size(), c.size());
        // ... calculate local_element_matrix values (e.g., all 1.0 for the test) ...
        local_element_matrix = Kratos::Matrix(c.size(), c.size(), 1.0); // Example from test

        A_glob_assem.Assemble(local_element_matrix, c);
    }
    A_glob_assem.FinalizeAssemble();
    // Check some values based on adding 1.0 matrices for {0,1} and {1,2}
    // A(0,0) should be 1.0, A(0,1) should be 1.0
    // A(1,0) should be 1.0, A(1,1) should be 1.0 (from {0,1}) + 1.0 (from {1,2}) = 2.0
    // A(1,2) should be 1.0, A(2,1) should be 1.0, A(2,2) should be 1.0
    KRATOS_EXPECT_EQ(A_glob_assem(1,1), 2.0); // Example check
    ```

*   `template<class TMatrixType, class TIndexVectorType> void Assemble(const TMatrixType& rMatrixInput, const TIndexVectorType& rRowEquationId, const TIndexVectorType& rColEquationId)`: Assembles a local (possibly rectangular) matrix `rMatrixInput` into the global `CsrMatrix`.
    *   `rRowEquationId`: Maps local rows of `rMatrixInput` to global rows of `*this`.
    *   `rColEquationId`: Maps local columns of `rMatrixInput` to global columns of `*this`.

    ```cpp
    // From KratosCoreFastSuite :: CSRMatrixRectangularMatrix (conceptual)
    // Assume CsrMatrix 'A_rect_assem' is constructed for a rectangular structure.
    Kratos::SparseContiguousRowGraph<> Agraph_rect_assem(2); // 2 rows
    // Define some sparsity, e.g. A(0,0), A(0,1), A(1,1), A(1,2)
    Agraph_rect_assem.AddEntry(0,0); Agraph_rect_assem.AddEntry(0,1);
    Agraph_rect_assem.AddEntry(1,1); Agraph_rect_assem.AddEntry(1,2);
    Agraph_rect_assem.Finalize();
    Kratos::CsrMatrix<double> A_rect_assem(Agraph_rect_assem);
    A_rect_assem.SetColSize(3); // 3 columns

    std::vector<Kratos::CsrMatrix<double>::IndexType> row_ids_rect = {0, 1}; // Global rows
    std::vector<Kratos::CsrMatrix<double>::IndexType> col_ids_rect = {1, 2}; // Global columns
    Kratos::Matrix local_contribution_matrix(row_ids_rect.size(), col_ids_rect.size()); // 2x2
    local_contribution_matrix(0,0) = 5.0; // Contributes to A(row_ids_rect[0], col_ids_rect[0]) = A(0,1)
    local_contribution_matrix(0,1) = 6.0; // Contributes to A(row_ids_rect[0], col_ids_rect[1]) = A(0,2)
    local_contribution_matrix(1,0) = 7.0; // Contributes to A(row_ids_rect[1], col_ids_rect[0]) = A(1,1)
    local_contribution_matrix(1,1) = 8.0; // Contributes to A(row_ids_rect[1], col_ids_rect[1]) = A(1,2)

    A_rect_assem.Assemble(local_contribution_matrix, row_ids_rect, col_ids_rect);
    KRATOS_EXPECT_EQ(A_rect_assem(0,1), 5.0); // Example check
    KRATOS_EXPECT_EQ(A_rect_assem(1,2), 8.0); // Example check
    ```

*   `void AssembleEntry(const TDataType Value, const IndexType GlobalI, const IndexType GlobalJ)`: Atomically adds `Value` to the entry at global row `GlobalI` and global column `GlobalJ`. The entry `(GlobalI, GlobalJ)` must already exist in the matrix's sparsity pattern.
    This is useful for adding individual contributions or for constructing the matrix if the pattern is known and entries are added one by one (though `Assemble` is often more efficient for FEM).

    ```cpp
    // From KratosCoreFastSuite :: CSRMatrixSmallRectangularMatrixMultiply (for matrix A)
    // Assume CsrMatrix 'A_entry_assem' is constructed from 'Agraph_entry_assem' defining its pattern.
    Kratos::SparseContiguousRowGraph<> Agraph_entry_assem(3);
    Agraph_entry_assem.AddEntry(0,0); Agraph_entry_assem.AddEntry(0,3); Agraph_entry_assem.AddEntry(0,4);
    Agraph_entry_assem.AddEntry(1,1);
    Agraph_entry_assem.AddEntry(2,3); Agraph_entry_assem.AddEntry(2,4);
    Agraph_entry_assem.Finalize();
    Kratos::CsrMatrix<double> A_entry_assem(Agraph_entry_assem);
    A_entry_assem.SetColSize(5); // Max column index is 4, so 5 columns

    A_entry_assem.BeginAssemble();
    A_entry_assem.AssembleEntry(1.0, 0, 0); // A_entry_assem(0,0) += 1.0
    A_entry_assem.AssembleEntry(7.0, 0, 3); // A_entry_assem(0,3) += 7.0
    A_entry_assem.AssembleEntry(2.0, 0, 4); // A_entry_assem(0,4) += 2.0
    A_entry_assem.AssembleEntry(3.0, 1, 1); // A_entry_assem(1,1) += 3.0
    A_entry_assem.AssembleEntry(7.0, 2, 3); // A_entry_assem(2,3) += 7.0
    A_entry_assem.AssembleEntry(7.0, 2, 4); // A_entry_assem(2,4) += 7.0
    A_entry_assem.FinalizeAssemble();
    KRATOS_EXPECT_EQ(A_entry_assem(0,3), 7.0); // Example check
    ```

### Other Operations

#### Boundary Conditions

*   `template<class TVectorType1, class TVectorType2=TVectorType1> void ApplyHomogeneousDirichlet(const TVectorType1& rFreeDofsVector, const TDataType DiagonalValue, TVectorType2& rRHS)`: Modifies the matrix `A` and a right-hand-side vector `rRHS` to apply homogeneous Dirichlet boundary conditions.
    *   `rFreeDofsVector`: A vector (e.g., `SystemVector`) of the same size as the matrix dimension. An entry `rFreeDofsVector[i] == 0.0` indicates that degree of freedom (DOF) `i` is fixed (Dirichlet). An entry `rFreeDofsVector[i] == 1.0` indicates it's a free DOF.
    *   `DiagonalValue`: The value to be set on the diagonal `A(i,i)` for rows `i` corresponding to fixed DOFs.
    *   `rRHS`: The right-hand side vector. `rRHS[i]` is set to `0.0` for fixed DOFs `i`.
    *   **Matrix Modification Details:**
        *   If DOF `i` is fixed (`rFreeDofsVector[i] == 0.0`):
            *   Row `i` of matrix `A` is zeroed out, i.e., `A(i,j) = 0.0` for all `j != i`.
            *   The diagonal element `A(i,i)` is set to `DiagonalValue`.
        *   If DOF `i` is free (`rFreeDofsVector[i] == 1.0`):
            *   For each entry `A(i,j)` in row `i`, it's multiplied by `rFreeDofsVector[j]`. This means if column `j` corresponds to a fixed DOF, `A(i,j)` becomes `0.0`, effectively removing the influence of fixed DOFs on free DOF equations.

    ```cpp
    // Conceptual Example:
    // Assume CsrMatrix 'A_bc' (e.g., 3x3) and SystemVector 'b_bc' are set up.
    Kratos::CsrMatrix<double>::MatrixMapType bc_map;
    bc_map[{0,0}]=2; bc_map[{0,1}]=-1;
    bc_map[{1,0}]=-1; bc_map[{1,1}]=2; bc_map[{1,2}]=-1;
    bc_map[{2,1}]=-1; bc_map[{2,2}]=2;
    Kratos::CsrMatrix<double> A_bc(bc_map);
    A_bc.ComputeColSize();

    Kratos::Vector b_bc(A_bc.size1());
    b_bc[0] = 1; b_bc[1] = 2; b_bc[2] = 3;
    // A_bc = [[2, -1, 0],
    //        [-1, 2, -1],
    //        [0, -1, 2]]  (Note: map doesn't add (0,2) or (2,0) if not specified)
    // b_bc = [1, 2, 3]

    Kratos::Vector free_dofs_bc(A_bc.size1());
    free_dofs_bc[0] = 1.0; // DOF 0 is free
    free_dofs_bc[1] = 0.0; // DOF 1 is fixed
    free_dofs_bc[2] = 1.0; // DOF 2 is free

    double diagonal_val_for_fixed_dofs = 1.0;
    A_bc.ApplyHomogeneousDirichlet(free_dofs_bc, diagonal_val_for_fixed_dofs, b_bc);

    // Expected result:
    // A_bc would become approximately:
    // [[2,  0,  0],  // Row 0: A(0,1) zeroed as DOF 1 is fixed
    //  [0,  1,  0],  // Row 1: Fixed DOF row, A(1,1) = diagonal_val_for_fixed_dofs, off-diagonals zeroed
    //  [0,  0,  2]]  // Row 2: A(2,1) zeroed as DOF 1 is fixed
    // b_bc would become:
    // [1, 0, 3]     // b_bc[1] (fixed DOF) is zeroed
    KRATOS_EXPECT_EQ(A_bc(1,1), diagonal_val_for_fixed_dofs); // Example check
    KRATOS_EXPECT_EQ(b_bc[1], 0.0); // Example check
    ```

#### Diagonal Properties

*   `TDataType NormDiagonal() const`: Calculates and returns the L2 norm (square root of sum of squares) of the diagonal elements of the matrix.
    ```cpp
    // From KratosCoreFastSuite :: CSRMatrixDiagonalValues
    // Assume CsrMatrix 'A_diag' is initialized.
    Kratos::CsrMatrix<double>::MatrixMapType diag_map;
    diag_map[{0,0}]=6.0; diag_map[{1,1}]=1.0; diag_map[{2,2}]=2.0; diag_map[{3,3}]=4.0;
    diag_map[{4,4}]=5.0; diag_map[{5,5}]=6.0; diag_map[{6,6}]=1.0; diag_map[{7,7}]=2.0;
    diag_map[{8,8}]=3.0; diag_map[{9,9}]=4.0; diag_map[{10,10}]=5.0; diag_map[{11,11}]=6.0;
    diag_map[{12,12}]=1.0; diag_map[{13,13}]=2.0; diag_map[{14,14}]=3.0; diag_map[{15,15}]=4.0;
    diag_map[{16,16}]=5.0; diag_map[{17,17}]=6.0; diag_map[{18,18}]=1.0; diag_map[{19,19}]=2.0;
    diag_map[{20,20}]=3.0; diag_map[{21,21}]=4.0; diag_map[{22,22}]=5.0; diag_map[{23,23}]=6.0;
    Kratos::CsrMatrix<double> A_diag(diag_map); // Using map from test data for example
    double norm_diag = A_diag.NormDiagonal();
    KRATOS_CHECK_NEAR(norm_diag, 22.135943621179, 1.0e-12); // Example from test
    ```

*   `TDataType MaxDiagonal() const`: Returns the maximum absolute value found on the diagonal of the matrix.
    ```cpp
    // From KratosCoreFastSuite :: CSRMatrixDiagonalValues
    // Assume CsrMatrix 'A_diag' (from above) is initialized.
    double max_diag = A_diag.MaxDiagonal();
    KRATOS_CHECK_NEAR(max_diag, 6.0, 1.0e-12); // Example from test
    ```

*   `TDataType MinDiagonal() const`: Returns the minimum absolute value found on the diagonal of the matrix.
    ```cpp
    // From KratosCoreFastSuite :: CSRMatrixDiagonalValues
    // Assume CsrMatrix 'A_diag' (from above) is initialized.
    double min_diag = A_diag.MinDiagonal();
    KRATOS_CHECK_NEAR(min_diag, 1.0, 1.0e-12); // Example from test
    ```

### Input and output

These methods provide ways to get textual information about the matrix or print its contents, which is primarily useful for debugging.

*   `std::string Info() const`: Returns a short string with basic information about the `CsrMatrix` class.
    ```cpp
    Kratos::CsrMatrix<double> A_info;
    std::string matrix_info = A_info.Info();
    std::cout << matrix_info << std::endl; 
    // Output would be something like: "CsrMatrix"
    ```

*   `void PrintInfo(std::ostream& rOStream) const`: Prints basic information about the matrix to the provided output stream. This typically includes the matrix type.
    ```cpp
    Kratos::CsrMatrix<double> A_pinfo;
    A_pinfo.PrintInfo(std::cout); 
    // Output to console:
    // CsrMatrix
    ```

*   `void PrintData(std::ostream& rOStream) const`: Prints detailed data of the matrix to the output stream. This includes its dimensions (`size1`, `size2`), number of non-zeros (`nnz`), and the content of its internal CSR arrays (`index1_data`, `index2_data`, `value_data`).
    ```cpp
    // Assume CsrMatrix 'A_pdata' is initialized, e.g. a small 2x2 matrix:
    A_pdata(0,0) = 1.0, A_pdata(0,1) = 2.0
    A_pdata(1,1) = 3.0
    // (Sparsity pattern must allow these entries)

    // Example setup:
    Kratos::CsrMatrix<double>::MatrixMapType data_map_pdata;
    data_map_pdata[{0,0}] = 1.0; data_map_pdata[{0,1}] = 2.0;
    data_map_pdata[{1,1}] = 3.0;
    Kratos::CsrMatrix<double> A_pdata(data_map_pdata);
    A_pdata.ComputeColSize(); // Ensure size2 is correct

    A_pdata.PrintData(std::cout);

    // Expected output (actual formatting might vary slightly):
    // size1 : 2
    // size2 : 2
    // nnz : 3
    // index1_data : 
    // 0,2,3, 
    // index2_data : 
    // 0,1,1, 
    // value_data  : 
    // 1,2,3, 
    ```

### Stream operators

*   `template<class TDataType> std::ostream& operator << (std::ostream& rOStream, const CsrMatrix<TDataType>& rThis)`:
    This overloaded output stream operator allows printing the matrix directly using `std::cout` or any other `std::ostream`. It typically calls `PrintInfo()` followed by `PrintData()`.

    ```cpp
    // Assume CsrMatrix 'A_print_stream' is initialized as in the PrintData example:
    Kratos::CsrMatrix<double>::MatrixMapType data_map_stream;
    data_map_stream[{0,0}] = 1.0; data_map_stream[{0,1}] = 2.0;
    data_map_stream[{1,1}] = 3.0;
    Kratos::CsrMatrix<double> A_print_stream(data_map_stream);
    A_print_stream.ComputeColSize(); // Ensure size2 is correct

    std::cout << "Printing matrix A_print_stream using stream operator:" << std::endl;
    std::cout << A_print_stream << std::endl;

    // Expected output (actual formatting might vary slightly):
    // Printing matrix A_print_stream using stream operator:
    // CsrMatrix
    //
    // size1 : 2
    // size2 : 2
    // nnz : 3
    // index1_data : 
    // 0,2,3, 
    // index2_data : 
    // 0,1,1, 
    // value_data  : 
    // 1,2,3, 
    ```

*   `template<class TDataType> std::istream& operator >> (std::istream& rIStream, CsrMatrix<TDataType>& rThis)`:
    The input stream operator is defined but is generally not used for populating a `CsrMatrix` due to the complexity of CSR structure. Constructing from a graph, map, or by reserving and assembling are the preferred methods.