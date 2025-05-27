---
title: DistributedCsrMatrix
keywords: DistributedCsrMatrix, Sparse Matrix, CSR, Distributed Computing, MPI, Parallel Linear Algebra
tags: [DistributedCsrMatrix, Data Structures, Linear Algebra, Distributed, MPI, HPC]
sidebar: kratos_for_developers
summary: Comprehensive documentation for the Kratos DistributedCsrMatrix class, a distributed Compressed Sparse Row (CSR) matrix designed for high-performance parallel computations using MPI.
---

## DistributedCsrMatrix Class

The `DistributedCsrMatrix` class in Kratos implements a distributed Compressed Sparse Row (CSR) matrix. This data structure is fundamental for parallel numerical simulations where large sparse matrices are partitioned and distributed across multiple processes, typically orchestrated via the Message Passing Interface (MPI). It extends the concepts of the serial `CsrMatrix` by incorporating mechanisms for managing data distribution, inter-process communication, and parallel assembly.

The `DistributedCsrMatrix` is crucial for large-scale simulations where a single process cannot hold the entire matrix in memory or when parallel processing is essential to achieve computational efficiency. It manages local sub-matrices (referred to as diagonal and off-diagonal blocks) and coordinates operations across these distributed parts.

For serial version see [CsrMatrix](CsrMatrix.md).

### Template Parameters

*   `TDataType`: The data type for the values stored in the matrix (e.g., `double`, `float`). Defaults to `double`.
*   `TIndexType`: The data type for row and column indices (e.g., `std::size_t`, `int`). Defaults to `std::size_t`.

### Type Definitions

*   `MpiIndexType`: Alias for `int`. This type is often used in MPI communication contexts, typically for counts or ranks.
*   `BlockMatrixType`: Alias for `CsrMatrix<TDataType, TIndexType>`. This is the type of the local matrices returned by `GetDiagonalBlock()` and `GetOffDiagonalBlock()`.
*   `MatrixMapType`: Alias for `std::map<std::pair<TIndexType, TIndexType>, TDataType>`. Used by the `ToMap()` method.

### Key Concepts in Distributed Environment

Understanding these concepts is vital for using `DistributedCsrMatrix` effectively:

*   **DataCommunicator**: An abstraction layer (typically wrapping MPI communicators like `MPI_Comm`) that manages communication (sending/receiving data) between different processes. Each `DistributedCsrMatrix` is associated with a `DataCommunicator`.
*   **RowNumbering / ColNumbering (`DistributedNumbering`)**: These objects define the mapping between global row/column indices and their local counterparts on each process. They also determine which process "owns" which global rows/columns, which is key for data locality and assembly.
*   **Diagonal Block**: The portion of the distributed matrix containing entries `(i,j)` where both the global row `i` and global column `j` are considered local to the current process, based on its `RowNumbering` and `ColNumbering`. This block is stored as a local `CsrMatrix`.
*   **Off-Diagonal Block**: The portion of the distributed matrix containing entries `(i,j)` where the global row `i` is local to the current process, but the global column `j` is *not* local. This block is also stored as a local `CsrMatrix`. Its column indices are local representations of these non-local global columns.
*   **DistributedVectorImporter / DistributedVectorExporter**: Utility classes crucial for SpMV-like operations.
    *   `DistributedVectorImporter` (`mpVectorImporter`): Gathers data from other processes. For `A*x`, if `A` has off-diagonal entries, the process needs the components of `x` corresponding to those non-local columns.
    *   `DistributedVectorExporter`: Scatters and sums data to other processes. For `A^T*x`, contributions to the result vector `y` might correspond to rows owned by other processes.
*   **Assembly**: The process of constructing the matrix by adding contributions (e.g., from element stiffness matrices in FEM).
    *   Entries `(GlobalI, GlobalJ)` where `GlobalI` is a local row are assembled directly into `DiagonalBlock` or `OffDiagonalBlock`.
    *   Contributions to entries where `GlobalI` is a non-local row are temporarily stored in a cache (`mNonLocalData`). These are sent to the owning process during `FinalizeAssemble()`.
*   **Communication Caches for Assembly**:
    *   `mNonLocalData`: A flat data array storing values to be sent to other processes.
    *   `mSendCachedIJ`, `mRecvCachedIJ`: Store pairs of global (row, col) indices for non-local entries that this process contributes to or receives contributions for.
    *   `mPointersToSendValues`, `mPointersToRecvValues`: Vectors of pointers providing direct access to locations within `mNonLocalData` (for sending) or the local block matrices (for receiving/assembling contributions from other processes). These are reconstructed after serialization.
    *   `mfem_assemble_colors`: Used for communication scheduling (`MPIColoringUtilities`) to manage non-blocking sends and receives efficiently during `FinalizeAssemble()`, avoiding deadlocks and improving performance.

### Lifecycle Methods

#### Constructors

*   `DistributedCsrMatrix()`: Default constructor. Creates an empty, uninitialized matrix. A `DataCommunicator` and graph structure must be provided later for the matrix to be functional.
    ```cpp
    Kratos::DistributedCsrMatrix<double> A;
    // A is not usable until properly initialized, e.g., via move assignment or by setting up its graph.
    ```

*   `DistributedCsrMatrix(const DataCommunicator& rComm)`: Constructor taking a `DataCommunicator`. The matrix is associated with this communicator, but its sparsity pattern remains undefined.
    ```cpp
    const Kratos::DataCommunicator& world_comm = Kratos::ParallelEnvironment::GetDefaultDataCommunicator();
    Kratos::DistributedCsrMatrix<double> A(world_comm);
    // A is associated with world_comm but has no entries or defined structure.
    ```

*   `DistributedCsrMatrix(const DistributedSparseGraph<TIndexType>& rSparseGraph)`: The primary constructor. It creates a `DistributedCsrMatrix` based on a pre-defined `DistributedSparseGraph`, which specifies the sparsity pattern and data distribution.
    *   The `DataCommunicator` is obtained from `rSparseGraph`.
    *   `RowNumbering` is copied from the graph.
    *   `ColNumbering` is established: if the matrix is globally square, column partitioning usually mirrors row partitioning. Otherwise, columns are distributed based on their total count.
    *   The `DiagonalBlock` and `OffDiagonalBlock` `CsrMatrix` structures are built by counting local and non-local entries for each local row.
    *   Internal mappings (`mOffDiagonalLocalIds`, `mOffDiagonalGlobalIds`) are created to translate between global column indices of non-local entries and their local representation in the `OffDiagonalBlock`.
    *   CSR arrays (`index1_data`, `index2_data`, `value_data`) for both blocks are allocated, column indices are sorted, and values are initialized to zero.
    *   Communication caches for assembly (`PrepareNonLocalCommunications`) and the `DistributedVectorImporter` (`mpVectorImporter`) for SpMV operations are prepared.

    ```cpp
    // Simplified example: Constructing a 5x5 DistributedCsrMatrix distributed over processes.
    KRATOS_TRY

    const Kratos::DataCommunicator& comm = Kratos::ParallelEnvironment::GetDefaultDataCommunicator();
    int rank = comm.Rank();
    int world_size = comm.Size();
    using IndexType = std::size_t;

    // 1. Define row partitioning (e.g., for a 5x5 matrix)
    // This defines how many rows each rank owns.
    std::vector<IndexType> row_partition_sizes(world_size);
    IndexType total_rows = 5;
    IndexType rows_per_rank_base = total_rows / world_size;
    IndexType remainder_rows = total_rows % world_size;
    for (int i = 0; i < world_size; ++i) {
        row_partition_sizes[i] = rows_per_rank_base + (i < remainder_rows ? 1 : 0);
    }
    auto p_row_numbering = Kratos::make_shared<Kratos::DistributedNumbering<IndexType>>(comm, row_partition_sizes);

    // 2. Create and populate a DistributedSparseGraph
    Kratos::DistributedSparseGraph<IndexType> dist_graph(*p_row_numbering);

    // Each rank adds entries ONLY for its local rows.
    // Example: A simple tridiagonal pattern
    for (IndexType i_global = p_row_numbering->MinId(); i_global <= p_row_numbering->MaxId(); ++i_global) {
        dist_graph.AddEntry(i_global, i_global); // Diagonal
        if (i_global > 0) {
            dist_graph.AddEntry(i_global, i_global - 1); // Lower band
        }
        if (i_global < total_rows - 1) {
            dist_graph.AddEntry(i_global, i_global + 1); // Upper band
        }
    }
    dist_graph.Finalize();

    // 3. Construct the DistributedCsrMatrix
    Kratos::DistributedCsrMatrix<double, IndexType> A(dist_graph);

    // Matrix A is now initialized with the sparsity of dist_graph, values are 0.0.
    A.GetDiagonalBlock() and A.GetOffDiagonalBlock() are set up.
    KRATOS_INFO("DistributedCsrMatrixExample") << "Rank " << rank << " constructed matrix A. Local rows: " << A.local_size1() << std::endl;

    KRATOS_CATCH("")
    ```

*   `explicit DistributedCsrMatrix(const DistributedCsrMatrix& rOtherMatrix)`: Copy constructor. Performs a deep copy of all components: `DataCommunicator` (pointer), numbering objects, block matrices (`DiagonalBlock`, `OffDiagonalBlock`), and all communication caches including the `DistributedVectorImporter`.
    ```cpp
    Kratos::DistributedCsrMatrix<double, IndexType> matrix_A(dist_graph); // from previous example
    Kratos::DistributedCsrMatrix<double, IndexType> matrix_B(matrix_A);
    // matrix_B is an independent, deep copy of matrix_A.
    ```

*   `DistributedCsrMatrix(DistributedCsrMatrix&& rOtherMatrix)`: Move constructor. Efficiently transfers ownership of resources from `rOtherMatrix` to the newly created matrix. `rOtherMatrix` is left in a valid but typically empty or undefined state.
    ```cpp
    Kratos::DistributedCsrMatrix<double, IndexType> matrix_C(dist_graph); // from previous example
    Kratos::DistributedCsrMatrix<double, IndexType> matrix_D(std::move(matrix_C));
    matrix_D now holds the data previously in matrix_C. matrix_C is moved-from.
    ```

#### Destructor

*   `~DistributedCsrMatrix()`: Destructor. Releases resources, including unique pointers to numbering objects, block matrices, and the vector importer.

#### Assignment Operators

*   `DistributedCsrMatrix& operator=(const DistributedCsrMatrix& rOther) = delete;`: Copy assignment is deleted to prevent accidental expensive copies. Use the copy constructor for explicit duplication.

*   `DistributedCsrMatrix& operator=(DistributedCsrMatrix&& rOtherMatrix)`: Move assignment operator. Transfers ownership of resources from `rOtherMatrix`.
    ```cpp
    Kratos::DistributedCsrMatrix<double, IndexType> matrix_E(dist_graph_E);
    Kratos::DistributedCsrMatrix<double, IndexType> matrix_F(dist_graph_F);
    matrix_E = std::move(matrix_F); // matrix_E now owns data from matrix_F.
    ```

### Main Operators and Member Functions

#### General Utility and Data Access

*   `void Clear()`: This method is present in the header but **its functionality to fully clear and reset a distributed matrix is not completely implemented or straightforward**. Clearing a distributed matrix would involve resetting local blocks and potentially complex re-initialization of communication structures. For a full reset, re-constructing from a graph is often safer.

*   `const DistributedNumbering<TIndexType>& GetRowNumbering() const`: Returns a const reference to the row `DistributedNumbering` object.
*   `const DistributedNumbering<TIndexType>& GetColNumbering() const`: Returns a const reference to the column `DistributedNumbering` object.
*   `typename DistributedNumbering<TIndexType>::UniquePointer& pGetRowNumbering()`: Returns a reference to the `std::unique_ptr` managing the row numbering.
*   `typename DistributedNumbering<TIndexType>::UniquePointer& pGetColNumbering()`: Returns a reference to the `std::unique_ptr` managing the column numbering.

*   `void SetValue(const TDataType Value)`: Sets all *existing* non-zero entries in both `DiagonalBlock` and `OffDiagonalBlock` to the specified `Value`. This does not change the sparsity pattern.
    ```cpp
    A.SetValue(0.0); // Zeroes out all structurally non-zero entries.
    ```

*   `TIndexType local_size1() const`: Returns the number of rows managed by the current process (i.e., `GetDiagonalBlock().size1()`).
*   `TIndexType size2() const`: Returns the total number of columns in the global matrix (`GetColNumbering().Size()`).
*   `TIndexType local_nnz() const`: Returns the number of stored entries in the `DiagonalBlock` only. To get total local stored entries: `A.GetDiagonalBlock().nnz() + A.GetOffDiagonalBlock().nnz()`.

*   `const DataCommunicator& GetComm() const`: Returns a const reference to the associated `DataCommunicator`.
*   `const DataCommunicator* pGetComm() const`: Returns a pointer to the `DataCommunicator`.

*   `BlockMatrixType& GetDiagonalBlock()`: Returns a writable reference to the local `CsrMatrix` storing the diagonal part.
*   `const BlockMatrixType& GetDiagonalBlock() const`: Const version.
*   `BlockMatrixType& GetOffDiagonalBlock()`: Returns a writable reference to the local `CsrMatrix` storing the off-diagonal part.
*   `const BlockMatrixType& GetOffDiagonalBlock() const`: Const version.
    ```cpp
    Kratos::CsrMatrix<double, IndexType>& local_diag_A = A.GetDiagonalBlock();
    Kratos::CsrMatrix<double, IndexType>& local_offdiag_A = A.GetOffDiagonalBlock();
    local_diag_A(local_row_idx, local_col_idx) = new_value; // Direct manipulation (use with care)
    ```

*   `const std::map<TIndexType, TIndexType>& GetOffDiagonalLocalIds() const`: Returns a map where keys are global column indices of non-local entries, and values are their corresponding local column indices within the `OffDiagonalBlock`.
*   `const DenseVector<TIndexType>& GetOffDiagonalGlobalIds() const`: Returns a vector where `mOffDiagonalGlobalIds[local_col_idx]` gives the global column index for `local_col_idx` in the `OffDiagonalBlock`.

*   `TIndexType GetOffDiagonalBlockLocalId(TIndexType GlobalJ) const`: Converts a non-local global column index `GlobalJ` to its local index in the `OffDiagonalBlock`. Throws an error if `GlobalJ` is not part of the off-diagonal pattern.
*   `TIndexType GetOffDiaGlobalId(TIndexType LocalJInOffDiag) const`: Converts a local column index `LocalJInOffDiag` from the `OffDiagonalBlock` back to its original global column index.

*   `TDataType& GetLocalDataByGlobalId(TIndexType GlobalI, TIndexType GlobalJ)`: Accesses an entry `(GlobalI, GlobalJ)` for modification. `GlobalI` *must* be a local row.
    *   If `GlobalJ` is a local column (i.e., `GetColNumbering().IsLocal(GlobalJ)` is true), it accesses `DiagonalBlock`.
    *   If `GlobalJ` is a non-local column, it accesses `OffDiagonalBlock` (using internal mapping).
    *   This is crucial for assembly of contributions to local rows.
    ```cpp
    // If GlobalI is a row owned by current process:
    if (A.GetRowNumbering().IsLocal(GlobalI)) {
        if (A.GetColNumbering().IsLocal(GlobalJ)) { // Column is also local
            A.GetLocalDataByGlobalId(GlobalI, GlobalJ) = 5.0;
        } else { // Column is non-local
            // This assumes (GlobalI, GlobalJ) is part of the sparsity pattern
            A.GetLocalDataByGlobalId(GlobalI, GlobalJ) = 6.0;
        }
    }
    ```

*   `TDataType& GetNonLocalCachedDataByGlobalId(TIndexType GlobalI, TIndexType GlobalJ)`: Accesses an entry `(GlobalI, GlobalJ)` in the temporary send cache (`mNonLocalData`). This is used during assembly for contributions where row `GlobalI` is *not* local to the current process. These cached values are communicated to the owner process during `FinalizeAssemble()`. `GlobalI` *must* be non-local.
    ```cpp
    // If GlobalI is a row NOT owned by current process:
    if (!A.GetRowNumbering().IsLocal(GlobalI)) {
        // This assumes (GlobalI, GlobalJ) has been prepared for non-local communication
        A.GetNonLocalCachedDataByGlobalId(GlobalI, GlobalJ) += contribution_to_non_local_row;
    }
    ```

*   `DenseVector<TIndexType> GetDiagonalIndex2DataInGlobalNumbering() const`: Returns a `DenseVector` containing the global column indices for entries in the `DiagonalBlock`, ordered as they appear in the block's `index2_data`.
*   `DenseVector<TIndexType> GetOffDiagonalIndex2DataInGlobalNumbering() const`: Similar to above, but for the `OffDiagonalBlock`, returning global column indices.

#### Assembly Operations

Assembly updates matrix entries with contributions, typically from element matrices in Finite Element Method (FEM) simulations. The sparsity pattern (graph) must be fixed before assembly begins. Atomic operations (`#pragma omp atomic`) are used internally for thread-safe updates to matrix values during assembly.

*   `void BeginAssemble()`: Initializes the assembly. Primarily, this zeroes out `mNonLocalData` (the cache for contributions to non-local rows).

*   `void FinalizeAssemble()`: Completes the assembly process. This is a collective operation involving communication:
    1.  It iterates through communication "colors" (from `mfem_assemble_colors`, determined by `MPIColoringUtilities::ComputeCommunicationScheduling`). Coloring optimizes non-blocking communication by grouping messages to avoid deadlocks and contention.
    2.  For each color:
        *   Data from `mNonLocalData` (contributions to non-local rows, accessed via `mPointersToSendValues`) is prepared for sending.
        *   `DataCommunicator::SendRecv` exchanges this data: each process sends its contributions and receives contributions from other processes for rows it owns.
        *   Received data (accessed via `mPointersToRecvValues`) is atomically added to the appropriate entries in the local `DiagonalBlock` or `OffDiagonalBlock`.

*   `template<class TMatrixType, class TIndexVectorType> void Assemble(const TMatrixType& rMatrixInput, const TIndexVectorType& rEquationId)`: Assembles a square local matrix `rMatrixInput` into the distributed matrix.
    *   `rMatrixInput`: The local matrix (e.g., element stiffness matrix).
    *   `rEquationId`: A vector mapping local indices `(0..N-1)` of `rMatrixInput` to global row AND column indices.
    *   For each `rMatrixInput(i_local, j_local)`:
        *   Global row `GlobalI = rEquationId[i_local]`, Global column `GlobalJ = rEquationId[j_local]`.
        *   If `GlobalI` is local: `rMatrixInput(i_local, j_local)` is atomically added to `GetLocalDataByGlobalId(GlobalI, GlobalJ)`.
        *   If `GlobalI` is non-local: `rMatrixInput(i_local, j_local)` is atomically added to `GetNonLocalCachedDataByGlobalId(GlobalI, GlobalJ)`.

*   `void AssembleEntry(const TDataType Value, const TIndexType GlobalI, const TIndexType GlobalJ)`: Assembles a single `Value` at global entry `(GlobalI, GlobalJ)`.
    *   If `GlobalI` is local: `Value` is atomically added to `GetLocalDataByGlobalId(GlobalI, GlobalJ)`.
    *   If `GlobalI` is non-local: `Value` is atomically added to `GetNonLocalCachedDataByGlobalId(GlobalI, GlobalJ)`.

*   `template<class TMatrixType, class TIndexVectorType> void Assemble(const TMatrixType& rMatrixInput, const TIndexVectorType& rRowEquationId, const TIndexVectorType& rColEquationId)`: Assembles a (possibly rectangular) local matrix `rMatrixInput`.
    *   `rRowEquationId`: Maps local rows of `rMatrixInput` to global rows.
    *   `rColEquationId`: Maps local columns of `rMatrixInput` to global columns.
    *   Logic is analogous to the square `Assemble`.

    ```cpp
    // Conceptual Assembly Workflow:
    // Assume 'A' is an initialized DistributedCsrMatrix.
    KRATOS_TRY

    A.BeginAssemble(); // Zero out non-local cache

    // --- Example: Assembling a single entry ---
    IndexType global_row_idx = ...; // Determine global row index
    IndexType global_col_idx = ...; // Determine global col index
    double value_to_add = 1.23;

    if (A.GetRowNumbering().IsLocal(global_row_idx)) {
        // Check if (global_row_idx, global_col_idx) is in the sparsity pattern of local blocks
        // This check is implicitly handled by GetLocalDataByGlobalId if the graph was built correctly.
        A.AssembleEntry(value_to_add, global_row_idx, global_col_idx);
    } else {
        // This assumes (global_row_idx, global_col_idx) was prepared for non-local communication
        // (i.e., it's in mSendCachedIJ for this rank).
        A.AssembleEntry(value_to_add, global_row_idx, global_col_idx);
    }

    // --- Example: Assembling a local dense matrix (e.g., from an element) ---
    Kratos::Matrix local_matrix(3,3); // Example 3x3 element matrix
    local_matrix(0,0) = 1.0; local_matrix(0,1) = 0.1; /* ... fill ... */
    std::vector<IndexType> equation_ids = {10, 12, 15}; // Global DOFs for this element

    A.Assemble(local_matrix, equation_ids);

    A.FinalizeAssemble(); // Communicate and sum all non-local contributions

    KRATOS_CATCH("")
    ```

#### Matrix-Vector Products (SpMV)

These operations perform `y = alpha*A*x + beta*y` or its transpose. They are fundamental in iterative solvers.

*   `void SpMV(const DistributedSystemVector<TDataType,TIndexType>& rX, DistributedSystemVector<TDataType,TIndexType>& rY) const`: Performs **rY += A\*rX**.
    *   `rX`: A `DistributedSystemVector` compatible with the matrix's column numbering.
    *   `rY`: A `DistributedSystemVector` compatible with the matrix's row numbering. `rY` is updated.
    *   **Operation Steps:**
        1.  `mpVectorImporter->ImportData(rX)`: Gathers non-local components of `rX` needed for the `OffDiagonalBlock` product. This involves MPI communication. The imported values are stored in `mOffDiagonalValuesCache`.
        2.  `GetOffDiagonalBlock().SpMV(mOffDiagonalValuesCache, rY.GetLocalData())`: Local SpMV of `OffDiagonalBlock` with imported `rX` data, adding to `rY`'s local data.
        3.  `GetDiagonalBlock().SpMV(rX.GetLocalData(), rY.GetLocalData())`: Local SpMV of `DiagonalBlock` with local `rX` data, adding to `rY`'s local data.

*   `void SpMV(const TDataType Alpha, const DistributedSystemVector<TDataType,TIndexType>& rX, const TDataType Beta, DistributedSystemVector<TDataType,TIndexType>& rY) const`: Performs **rY = Alpha\*A\*rX + Beta\*rY**.
    *   Incorporates scaling factors `Alpha` and `Beta`. If `Beta` is 0.0, `rY` is effectively overwritten. If `Alpha` is 1.0 and `Beta` is 1.0 (and `rY` was initially zero), it's `rY = A*rX`.
    *   The implementation scales `rY` by `Beta` first, then adds `Alpha * A * rX`.

    ```cpp
    // Assume 'A' is an assembled DistributedCsrMatrix.
    // Assume 'x_dist' and 'y_dist' are DistributedSystemVectors.
    KRATOS_TRY

    Kratos::DistributedSystemVector<double, IndexType> x_dist(A.GetColNumbering());
    Kratos::DistributedSystemVector<double, IndexType> y_dist(A.GetRowNumbering());

    // Initialize x_dist (e.g., with ones)
    for(IndexType i=0; i<x_dist.GetLocalData().size(); ++i) x_dist.GetLocalData()[i] = 1.0;
    x_dist.GetComm().Barrier(); // Ensure all ranks initialized before SpMV

    // Initialize y_dist (e.g., with zeros)
    for(IndexType i=0; i<y_dist.GetLocalData().size(); ++i) y_dist.GetLocalData()[i] = 0.0;
    y_dist.GetComm().Barrier();

    // Perform y_dist = 1.0 * A * x_dist + 0.0 * y_dist  (i.e., y_dist = A*x_dist)
    A.SpMV(1.0, x_dist, 0.0, y_dist);

    // y_dist now contains the result of A*x_dist on each process.
    KRATOS_INFO("SpMVExample") << "Rank " << A.GetComm().Rank() << " y_dist[0] = " << (y_dist.GetLocalData().size() > 0 ? y_dist.GetLocalData()[0] : 0.0) << std::endl;

    KRATOS_CATCH("")
    ```

*   `DistributedVectorExporter<TIndexType>* TransposeSpMV(const DistributedSystemVector<TDataType,TIndexType>& rX, DistributedSystemVector<TDataType,TIndexType>& rY, DistributedVectorExporter<TIndexType>* pTransposeExporter = nullptr) const`: Performs **rY += A<sup>T</sup>\*rX**.
    *   Transpose SpMV is more complex due to data dependencies.
    *   `GetDiagonalBlock().TransposeSpMV(rX.GetLocalData(), rY.GetLocalData())`: Local transpose SpMV for the diagonal block.
    *   `GetOffDiagonalBlock().TransposeSpMV(rX.GetLocalData(), mNonLocalTransposeValuesCache)`: Local transpose SpMV for the off-diagonal block. The results in `mNonLocalTransposeValuesCache` are contributions for `rY` components that are non-local (columns of original `A` become rows of `A^T`).
    *   `DistributedVectorExporter` (`pTransposeExporter` or a new one if `nullptr`): Manages summing these `mNonLocalTransposeValuesCache` contributions to the correct owning processes of `rY`. Creating an exporter can involve communication to set up patterns. Reusing an existing one is more efficient for repeated calls.
    *   Returns the exporter used (caller might need to manage its lifecycle if it was dynamically created by not passing one in).

*   `DistributedVectorExporter<TIndexType>* TransposeSpMV(TDataType Alpha, const DistributedSystemVector<TDataType,TIndexType>& rX, TDataType Beta, DistributedSystemVector<TDataType,TIndexType>& rY, DistributedVectorExporter<TIndexType>* pTransposeExporter = nullptr) const`: Performs **rY = Alpha\*A<sup>T</sup>\*rX + Beta\*rY**. Similar logic with scaling factors.

#### Norms and Diagonal Properties

These methods compute global properties of the matrix, involving communication.

*   `TDataType NormFrobenius() const`: Calculates the Frobenius norm: `sqrt(sum(A_ij^2))`. Each process computes the sum of squares for its local entries (diagonal and off-diagonal blocks). These local sums are then summed across all processes using `GetComm().SumAll()`, and the square root is taken.

*   `TDataType NormDiagonal() const`: Calculates the L2 norm of the main diagonal of the global matrix. Each process computes the sum of squares of its local diagonal entries (from `GetDiagonalBlock()`). These are summed via `GetComm().SumAll()`, and the square root is taken.

*   `TDataType MaxDiagonal() const`: Finds the maximum absolute value on the global main diagonal. Each process finds its local maximum diagonal entry, and `GetComm().MaxAll()` determines the global maximum.
*   `TDataType MinDiagonal() const`: Finds the minimum absolute value on the global main diagonal. Each process finds its local minimum, and `GetComm().MinAll()` determines the global minimum.

#### Other Methods

*   `MatrixMapType ToMap() const`: Converts the local part of the distributed matrix (both `DiagonalBlock` and `OffDiagonalBlock`) into a `std::map<std::pair<TIndexType, TIndexType>, TDataType>`, where keys are global `(row, column)` indices. This is a local operation; it does not gather the entire matrix.
    ```cpp
    Kratos::DistributedCsrMatrix<double, IndexType>::MatrixMapType local_map = A.ToMap();
    for(const auto& entry : local_map) {
        entry.first.first is global row, entry.first.second is global col, entry.second is value
    }
    ```

*   `typename CsrMatrix<TDataType,TIndexType>::Pointer ToSerialCSR(MpiIndexType TargetRank = 0) const`: Gathers all distributed parts of the matrix to a single process (`TargetRank`) and reconstructs it as a standard (serial) `CsrMatrix`.
    *   Each process serializes its `DiagonalBlock` and `OffDiagonalBlock` entries (global row, global col, value) into flat lists.
    *   `GetComm().Gatherv()` collects these lists on `TargetRank`.
    *   `TargetRank` then:
        1.  Builds a global `SparseContiguousRowGraph`.
        2.  Populates this graph with all received entries.
        3.  Creates a new `CsrMatrix` from this global graph and assembles the values.
    *   Returns a `shared_ptr` to the new serial `CsrMatrix` on `TargetRank`; other ranks return `nullptr`.
    *   **Caution**: This can be very memory-intensive on `TargetRank` for large matrices. Useful for debugging or interfacing with serial libraries.

    ```cpp
    Kratos::CsrMatrix<double, IndexType>::Pointer p_serial_A = A.ToSerialCSR(0); // Gather on rank 0
    if (A.GetComm().Rank() == 0 && p_serial_A) {
        // Rank 0 now has the full matrix in p_serial_A.
        std::cout << "Full matrix on Rank 0: " << (*p_serial_A) << std::endl;
    }
    ```

*   `void ApplyHomogeneousDirichlet(const DistributedSystemVector<TDataType, TIndexType>& rDofStatus, const TDataType rPrescribedValue, DistributedSystemVector<TDataType, TIndexType>& rRHS) const`: Modifies the matrix and a right-hand-side vector `rRHS` to apply homogeneous Dirichlet boundary conditions (typically setting diagonal entries to 1, off-diagonals in the row to 0, and adjusting RHS).
    *   Operates on `DiagonalBlock` similar to serial `CsrMatrix::ApplyDirichlet`.
    *   For `OffDiagonalBlock`, it uses `mpVectorImporter` to get non-local `rDofStatus` to correctly modify rows connected to fixed non-local DOFs.

### Input and Output

*   `std::string Info() const`: Returns a short string: `"DistributedCsrMatrix"`.
*   `void PrintInfo(std::ostream& rOStream) const`: Prints `"DistributedCsrMatrix"` to the given stream.
*   `void PrintData(std::ostream& rOStream) const`: Prints detailed information about the local data on the current process to `rOStream`. This includes:
    *   `DiagonalBlock`: sizes, nnz, index1_data, index2_data (local and global numbering), value_data.
    *   `OffDiagonalBlock`: sizes, nnz, index1_data, index2_data (local and global numbering), value_data.
    *   This method only shows data local to the calling rank. For the global matrix, use `ToSerialCSR()` and print the result.

### Serialization

Serialization allows saving and loading the `DistributedCsrMatrix` state.

*   `void save(Serializer& rSerializer) const`: Saves the matrix. Serialized components include:
    *   The name of the `DataCommunicator`.
    *   `RowNumbering` and `ColNumbering`.
    *   `DiagonalBlock` and `OffDiagonalBlock` (delegating to their own `save` methods).
    *   `mNonLocalData` (cache for non-local assembly contributions).
    *   Communication metadata: `mSendCachedIJ`, `mRecvCachedIJ`, `mfem_assemble_colors`.
    *   Mappings: `mOffDiagonalLocalIds`, `mOffDiagonalGlobalIds`.
    *   The `DistributedVectorImporter` (`mpVectorImporter`).
    *   **Note**: Raw pointer caches like `mPointersToSendValues` and `mPointersToRecvValues` are *not* saved directly.

*   `void load(Serializer& rSerializer)`: Loads the matrix. It deserializes all components saved by `save` and then calls `ReconstructDirectAccessVectors()` to rebuild the internal pointer-based caches (`mPointersToSendValues`, `mPointersToRecvValues`) essential for the assembly communication logic. This ensures the matrix is fully functional after loading.

### Example Usage (Conceptual High-Level Workflow)

```cpp
    // This is a high-level conceptual example of building and using a DistributedCsrMatrix.
    // Assume 'comm' is the DataCommunicator (e.g., Kratos::ParallelEnvironment::GetDefaultDataCommunicator()).
    // Assume 'model_part' is a Kratos ModelPart containing elements and conditions.
    KRATOS_TRY

    using IndexType = std::size_t;
    using DataType = double;

    const Kratos::DataCommunicator& comm = Kratos::ParallelEnvironment::GetDefaultDataCommunicator();

    // 1. Determine DOF distribution and create DistributedNumbering.
    //    This is typically done by a builder_and_solver or a utility.
    //    For this example, let's assume p_row_numbering is already created,
    //    defining how many DOFs (rows) each rank will manage.
    //    Kratos::shared_ptr<Kratos::DistributedNumbering<IndexType>> p_row_numbering = ...;
    //    For a concrete example of p_row_numbering, see the constructor example.
    //    If p_row_numbering is not available, this example cannot proceed.
    //    For placeholder:
    std::vector<IndexType> row_partition_sizes(comm.Size());
    IndexType total_rows_placeholder = 100; // Example total number of DOFs
    IndexType rows_per_rank_base_ph = total_rows_placeholder / comm.Size();
    IndexType remainder_rows_ph = total_rows_placeholder % comm.Size();
    for (int i = 0; i < comm.Size(); ++i) {
        row_partition_sizes[i] = rows_per_rank_base_ph + (i < remainder_rows_ph ? 1 : 0);
    }
    auto p_row_numbering = Kratos::make_shared<Kratos::DistributedNumbering<IndexType>>(comm, row_partition_sizes);

    // 2. Create DistributedSparseGraph: Iterate over local elements/conditions
    //    to define the matrix sparsity pattern.
    Kratos::DistributedSparseGraph<IndexType> graph(*p_row_numbering);

    // In a real scenario, you'd loop through elements/conditions:
    for (auto& elem : model_part.GetCommunicator().LocalMesh().Elements()) {
        std::vector<IndexType> elem_dof_ids;
        elem.GetDofList().GetEquationIds(elem_dof_ids); // Get global DOF IDs
        graph.AddEntries(elem_dof_ids); // Add connections for these DOFs
    }
    // For placeholder:
    for (IndexType i_global = p_row_numbering->MinId(); i_global <= p_row_numbering->MaxId(); ++i_global) {
        graph.AddEntry(i_global, i_global); // Diagonal entry
    }
    graph.Finalize();

    // 3. Create DistributedCsrMatrix.
    Kratos::DistributedCsrMatrix<DataType, IndexType> A(graph);

    // 4. Assemble the matrix.
    A.BeginAssemble();

   // In a real scenario, loop through elements/conditions again:
   Kratos::Matrix LHS_contribution; // Local element matrix
   std::vector<IndexType> GIDs;    // Global IDs for element DOFs
   for (auto& elem : model_part.GetCommunicator().LocalMesh().Elements()) {
       elem.CalculateLocalSystem(LHS_contribution, some_RHS_vector, process_info);
       elem.GetDofList().GetEquationIds(GIDs);
       A.Assemble(LHS_contribution, GIDs);
   }
    // For placeholder (assembling a diagonal value of 2.0 for each local row):
    for (IndexType i_global = p_row_numbering->MinId(); i_global <= p_row_numbering->MaxId(); ++i_global) {
        A.AssembleEntry(2.0, i_global, i_global);
    }
    A.FinalizeAssemble();

    // 5. Create DistributedSystemVectors for solution (x) and RHS (b).
    Kratos::DistributedSystemVector<DataType, IndexType> x(A.GetColNumbering());
    Kratos::DistributedSystemVector<DataType, IndexType> b(A.GetRowNumbering());

    //    Fill 'b' with values (e.g., from loads or body forces).
    //    Apply boundary conditions to A and b.
    //    (Example: x.Clear(); b.Clear(); /* ... fill b ... */ A.ApplyHomogeneousDirichlet(...); )

    // 6. Solve A*x = b using a distributed iterative solver.
    //    (e.g., Kratos::AMGCLSolver, Kratos::GMRESSolver with a distributed preconditioner)
    SomeSolverType solver;
    solver.Solve(A, x, b);

    // 7. 'x' now contains the distributed solution. Each process holds its local part of 'x'.
    const auto& local_solution_data = x.GetLocalData();
    KRATOS_INFO("DistributedWorkflow") << "Rank " << comm.Rank() << " solution (first local entry if any): "
                                        << (local_solution_data.size() > 0 ? local_solution_data[0] : 0.0) << std::endl;

    KRATOS_CATCH("")
```