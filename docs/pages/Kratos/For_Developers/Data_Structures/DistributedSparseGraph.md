---
title: DistributedSparseGraph
keywords: DistributedSparseGraph, Sparse Matrix, CSR, Distributed Computing, MPI, Parallel Linear Algebra
tags: [DistributedSparseGraph, Data Structures, Linear Algebra, Distributed, MPI, HPC]
sidebar: kratos_for_developers
summary: Comprehensive documentation for the Kratos DistributedSparseGraph class, a distributed sparse graph designed for high-performance parallel computations using MPI.
---

## DistributedSparseGraph Class

The `DistributedSparseGraph` class implements a distributed sparse graph data structure designed for parallel computations using MPI. It allows different MPI ranks to manage portions of a global graph, facilitating distributed assembly and storage of graph connectivity. This class is crucial for large-scale simulations where a single machine's memory is insufficient to hold the entire graph.

The graph is built upon a `DataCommunicator` which defines the MPI context. It internally creates and manages a `DistributedNumbering` object to handle the mapping between global and local row indices across ranks. Each rank stores its local portion of the graph (for rows it owns) and uses temporary storage for entries that belong to non-local rows, which are communicated during `Finalize()`.

**Important Note on Thread Safety:** While `DistributedSparseGraph` coordinates data across MPI ranks, operations that modify internal graph structures (especially those involving non-local entries before `Finalize()`) use locks for thread safety. However, the underlying local graph (`SparseContiguousRowGraph`) and non-local graph buffers (`SparseGraph`) might have their own considerations if populated directly in a multi-threaded fashion outside of the provided `AddEntry`/`AddEntries` methods. The header notes: "it is BY DESIGN NOT threadsafe! (a graph should be computed in each thread and then merged)" - this generally applies to how graph data might be assembled before being passed to this class or if manipulating its internal structures directly.

## Template Parameters

The `DistributedSparseGraph` class has one template parameter:

- `TIndexType`: Specifies the type of indices used for the graph nodes (e.g., `std::size_t`). This type should be consistent with the `DistributedNumbering` used.

## Lifecycle Methods

### Constructor

**`DistributedSparseGraph(DistributedNumbering<TIndexType>& rRowNumbering)`**
Constructs a `DistributedSparseGraph`.
- `rRowNumbering`: A reference to a `DistributedNumbering` object. This object defines how global row indices are distributed among MPI ranks and provides methods for converting between global and local indices. The `DataCommunicator` is obtained from this numbering.

```cpp
// Assuming r_comm is an MPI_Comm (e.g., MPI_COMM_WORLD)
Kratos::DataCommunicator& kratos_comm = Kratos::ParallelEnvironment::GetDataCommunicator("World");

// Example: Create a DistributedNumbering for rows
// Suppose total 100 rows, rank 0 owns 0-49, rank 1 owns 50-99 (simplified)
std::vector<Kratos::IndexType> partition;
// ... fill partition based on rank ...
// For rank 0: partition = {50} (meaning it owns 50 rows, indices 0 to 49 if ghost layers are 0)
// For rank 1: partition = {50} (meaning it owns 50 rows, indices 50 to 99 if ghost layers are 0)
// Actual partitioning depends on how DistributedNumbering is constructed.
// Let's assume a simple case for rRowNumbering construction:
Kratos::DistributedNumbering<Kratos::IndexType> row_numbering(kratos_comm, partition);

// Construct the DistributedSparseGraph
Kratos::DistributedSparseGraph<Kratos::IndexType> dist_graph(row_numbering);
```
### Destructor
The destructor handles the deallocation of resources managed by the `DistributedSparseGraph`.
```cpp
// dist_graph goes out of scope, its destructor is automatically called.
```
**Note:** `DistributedSparseGraph` does not explicitly define a copy constructor or assignment operator, implying default behavior or that they are not intended for typical use.

## Main Operators and Member Functions

### `AddEntry(const IndexType GlobalRowIndex, const IndexType GlobalColIndex)`
Adds a single directed edge (entry) from global row `GlobalRowIndex` to global column `GlobalColIndex`.
If `GlobalRowIndex` is owned by the current rank, the entry `(GlobalRowIndex, GlobalColIndex)` is added to the local part of the graph.
If `GlobalRowIndex` is owned by another rank, this information is stored temporarily and communicated during `Finalize()`.
```cpp
// Assuming dist_graph is initialized and rRowNumbering is available
// Kratos::IndexType my_rank = kratos_comm.Rank();
// Kratos::IndexType global_row_idx = /* some global row index */;
// Kratos::IndexType global_col_idx = /* some global column index */;

// if (row_numbering.IsLocal(global_row_idx)) { // Example: rank 0 owns global rows 0-4
//     dist_graph.AddEntry(global_row_idx, global_col_idx + my_rank); // Example entry
// }
// Example from test_distributed_sparse_graph.cpp:
// rank 0 adds entry (0,0) and (0,1)
// rank 1 adds entry (2,2) and (2,3)
// dist_graph.AddEntry(0,0); // If rank 0
// dist_graph.AddEntry(2,2); // If rank 1
```

### `AddEntries(const IndexType GlobalRowIndex, const std::vector<IndexType>& rGlobalColIndices)`
Adds multiple directed edges from a source global row `GlobalRowIndex` to all global column indices specified in `rGlobalColIndices`.
Similar to `AddEntry`, handles local and non-local `GlobalRowIndex` appropriately.
```cpp
// Kratos::IndexType global_row_idx = /* ... */;
// std::vector<Kratos::IndexType> global_target_cols = { /* ... */ };
// dist_graph.AddEntries(global_row_idx, global_target_cols);
// Example from test_distributed_sparse_graph.cpp:
// std::vector<int> indices_to_add = {0,1}; // connectivities for a specific element
// dist_graph.AddEntries(indices_to_add[0], indices_to_add); // Adds (0,0) and (0,1) if indices_to_add[0] is 0
```

### `AddEntries(const std::vector<IndexType>& rGlobalRowIndices, const std::vector<IndexType>& rGlobalColIndices)`
For each global row index `GR_i` in `rGlobalRowIndices`, adds directed edges from `GR_i` to all global column indices in `rGlobalColIndices`.
```cpp
// std::vector<Kratos::IndexType> global_source_rows = { /* ... */ };
// std::vector<Kratos::IndexType> global_target_cols = { /* ... */ };
// dist_graph.AddEntries(global_source_rows, global_target_cols);
```

### `AddEntries(const std::vector<std::vector<IndexType>>& rGlobalConnectivities)`
Adds entries based on a list of connectivities. Each inner vector in `rGlobalConnectivities` represents a clique of interconnected global indices (e.g., nodes of a finite element). For each index `I` in an inner vector, entries `(I,J)` are added for all other indices `J` in the same inner vector.
```cpp
// std::vector<std::vector<Kratos::IndexType>> connectivities;
// connectivities.push_back({0,1,2}); // Element 1 connecting nodes 0,1,2
// connectivities.push_back({2,3,0}); // Element 2 connecting nodes 2,3,0
// dist_graph.AddEntries(connectivities);
// For {0,1,2}: adds (0,0),(0,1),(0,2), (1,0),(1,1),(1,2), (2,0),(2,1),(2,2)
// (behavior regarding local/non-local rows applies)
```

### `Finalize()`
Finalizes the graph construction. This is a collective operation that must be called on all MPI ranks.
It performs the following key steps:
1. Exchanges information about non-local entries: If a rank tried to add an entry for a row it does not own, this information is sent to the owner rank.
2. Populates the `mNonLocalGraph` on each rank: This graph stores, for rows owned by the current rank, which other ranks also have these rows as non-local (ghost) rows and need to know about their connectivity.
```cpp
dist_graph.Finalize(); // Collective call
```

### `GetLocalGraph() const`
Returns a constant reference to the local part of the graph. This graph contains entries `(GlobalRowIndex, GlobalColIndex)` where `GlobalRowIndex` is owned by the current MPI rank.
The returned type is `const SparseGraph<TIndexType>&`.
```cpp
const Kratos::SparseGraph<Kratos::IndexType>& local_graph = dist_graph.GetLocalGraph();
// Iterate over local_graph to see locally owned rows and their connections
// for(auto it_row = local_graph.begin(); it_row != local_graph.end(); ++it_row) {
//     Kratos::IndexType global_row_idx = it_row.GetRowIndex();
//     const auto& neighbors = *it_row;
//     // ...
// }
```

### `GetNonLocalGraph() const`
Returns a constant reference to the non-local graph structure.
This is a `std::map<IndexType, std::vector<int>>`. The map key is a `GlobalRowIndex` (owned by the current rank). The associated `std::vector<int>` contains the ranks of the MPI processes that have this `GlobalRowIndex` as a ghost/non-local row and therefore need to be informed about its connectivity.
This is primarily for internal use or advanced scenarios.
```cpp
const std::map<Kratos::IndexType, std::vector<int>>& non_local_graph_info = dist_graph.GetNonLocalGraph();
// For a locally owned GlobalRowIndex, non_local_graph_info[GlobalRowIndex] lists ranks that ghost it.
```

### `GetRowNumbering() const`
Returns a constant reference to the `DistributedNumbering<TIndexType>` object used by the graph.
```cpp
const Kratos::DistributedNumbering<Kratos::IndexType>& numbering = dist_graph.GetRowNumbering();
// Kratos::IndexType owner_rank = numbering.OwnerRank(global_idx);
```

### `Size()`
Returns the total number of unique global row indices that have at least one entry in the graph, across all MPI ranks. This implies the "height" or number of active rows in the conceptual global matrix graph.
```cpp
Kratos::IndexType total_active_rows = dist_graph.Size();
```

### `LocalSize()`
Returns the number of unique global row indices owned by the current MPI rank that have at least one entry in the graph.
```cpp
Kratos::IndexType local_active_rows = dist_graph.LocalSize();
// This is equivalent to dist_graph.GetLocalGraph().Size() if GetLocalGraph().Size() counts non-empty rows.
// More precisely, it's the number of keys in the local graph part.
```

### `Has(const IndexType GlobalRowIndex, const IndexType GlobalColIndex) const`
Checks if a specific edge `(GlobalRowIndex, GlobalColIndex)` exists in the graph.
This operation is local if `GlobalRowIndex` is owned by the current rank. If `GlobalRowIndex` is owned by another rank, it currently might only check the sender's knowledge before `Finalize`, or it might not be well-defined for remote rows post-`Finalize` without further communication (the header suggests it checks `mpLocalGraph`). It's safest to assume this is primarily for locally owned rows.
```cpp
// Kratos::IndexType global_row_idx = /* ... */;
// Kratos::IndexType global_col_idx = /* ... */;
// if (dist_graph.GetRowNumbering().IsLocal(global_row_idx)) {
//    bool exists = dist_graph.Has(global_row_idx, global_col_idx);
//    // ...
// }
```

### `Clear()`
Clears all data in the `DistributedSparseGraph` on all ranks. This includes local entries, non-local communication buffers, and resets the finalized state.
This is a collective operation.
```cpp
dist_graph.Clear(); // Collective call
// dist_graph is now empty and ready to be rebuilt.
```

## Input and Output Methods

### `Info() const`
Returns a string containing basic information about the `DistributedSparseGraph`.
```cpp
std::cout << dist_graph.Info() << std::endl;
// Example output might include information about local/global size.
```

### `PrintInfo(std::ostream& rOStream) const`
Prints basic information about the `DistributedSparseGraph` to the specified output stream `rOStream`.
This is a collective operation and typically prints information only on the root rank (rank 0).
```cpp
// Only rank 0 prints:
if (dist_graph.GetRowNumbering().GetDataCommunicator().Rank() == 0) {
    dist_graph.PrintInfo(std::cout);
    std::cout << std::endl;
}
```

### `PrintData(std::ostream& rOStream) const`
Prints the data of the `DistributedSparseGraph` to the specified output stream `rOStream`.
This is a collective operation. Each rank prints its local data.
```cpp
// Each rank prints its local portion:
// std::cout << "Data on Rank " << dist_graph.GetRowNumbering().GetDataCommunicator().Rank() << ":" << std::endl;
// dist_graph.PrintData(std::cout);
// std::cout << std::endl;

// For controlled printing (e.g. rank by rank):
// for (int i = 0; i < kratos_comm.Size(); ++i) {
//     kratos_comm.Barrier();
//     if (kratos_comm.Rank() == i) {
//         std::cout << "Data for rank " << i << ":" << std::endl;
//         dist_graph.PrintData(std::cout);
//         std::cout << std::endl;
//     }
//     kratos_comm.Barrier();
// }
```

### Stream Operators
**`operator<<(std::ostream& rOStream, const DistributedSparseGraph<TIndexType>& rThis)`**
Overloads the `<<` operator for output streams. It typically calls `PrintInfo` and `PrintData`.
The behavior regarding which ranks print can be similar to `PrintInfo` or `PrintData`.
```cpp
// if (dist_graph.GetRowNumbering().GetDataCommunicator().Rank() == 0) {
//    std::cout << dist_graph << std::endl;
// }
```

## Serialization

The `DistributedSparseGraph` class provides `save` and `load` methods for serialization, allowing its state to be saved and restored. This is essential for checkpointing and restarting large simulations.

### `save(Serializer& rSerializer) const`
Serializes the `DistributedSparseGraph` object. This is a collective operation. Each rank serializes its local data.
```cpp
// Assuming rSerializer is a Kratos::Serializer object configured for distributed saving.
// dist_graph.save(rSerializer);
```

### `load(Serializer& rSerializer)`
Deserializes the `DistributedSparseGraph` object. This is a collective operation. Each rank deserializes its local data. The `DistributedNumbering` must be set up compatibly before loading.
```cpp
// Kratos::DistributedNumbering<Kratos::IndexType> row_numbering_for_load(...);
// Kratos::DistributedSparseGraph<Kratos::IndexType> graph_to_load(row_numbering_for_load);
// graph_to_load.load(rSerializer);
```

This documentation provides a foundational understanding of the `DistributedSparseGraph`. For specific details on advanced usage or the internals of communication patterns, direct inspection of the source code and test examples is recommended.
```md
# DistributedSparseGraph

## Overview

The `DistributedSparseGraph` class implements a distributed sparse graph data structure designed for parallel computations using MPI. It allows different MPI ranks to manage portions of a global graph, facilitating distributed assembly and storage of graph connectivity. This class is crucial for large-scale simulations where a single machine's memory is insufficient to hold the entire graph.

The graph is built upon a `DataCommunicator` which defines the MPI context and a `DistributedNumbering` which manages the mapping between global and local indices across ranks. Each rank stores its local portion of the graph and information about non-local entries (connections to rows owned by other ranks).

## Template Parameters

The `DistributedSparseGraph` class has one template parameter:

- `TIndexType`: Specifies the type of indices used for the graph nodes (e.g., `std::size_t`). This type should be consistent with the `DistributedNumbering` used.

## Lifecycle Methods

### Constructor

**`DistributedSparseGraph(const IndexType LocalSize, const DataCommunicator& rComm)`**
Constructs a `DistributedSparseGraph`.
- `LocalSize`: The number of rows that the current MPI rank will own.
- `rComm`: A reference to the `DataCommunicator` that defines the MPI context for this distributed graph.

The constructor initializes the local graph part (`mLocalGraph`) with `LocalSize`. It also internally creates a `DistributedNumbering` object (`mpRowNumbering`) based on `rComm` and the `LocalSize` provided by each rank. This internal numbering manages the global row indices and their distribution.

```cpp
const Kratos::DataCommunicator& kratos_comm = Kratos::ParallelEnvironment::GetDataCommunicator("World");
Kratos::IndexType num_rows_for_current_rank = 10; // Example: this rank will own 10 rows.
                                                // In practice, this value can differ per rank.

// Each rank constructs the graph specifying how many rows it owns.
// The DistributedNumbering is created internally to map these local sizes to a global indexing.
Kratos::DistributedSparseGraph<Kratos::IndexType> dist_graph(num_rows_for_current_rank, kratos_comm);
```
Example from `test_distributed_sparse_graph.cpp`:
```cpp
// dofs_bounds are computed per rank: dofs_bounds[0] is start_id, dofs_bounds[1] is end_id for this rank.
// So, dofs_bounds[1] - dofs_bounds[0] is the number of rows this rank owns.
// DistributedSparseGraph<IndexType> Agraph(dofs_bounds[1]-dofs_bounds[0], r_comm);
```

### Destructor
The destructor handles the deallocation of resources, including the internally created `DistributedNumbering`.
```cpp
// dist_graph goes out of scope, its destructor is automatically called.
```
**Copy and Assignment:**
The copy constructor and assignment operator for `DistributedSparseGraph` are explicitly deleted. This means instances of this class cannot be copied or assigned.
```cpp
// Kratos::DistributedSparseGraph<Kratos::IndexType> graph_copy(dist_graph); // Error: Copy constructor deleted
// Kratos::DistributedSparseGraph<Kratos::IndexType> another_graph(10, kratos_comm);
// another_graph = dist_graph; // Error: Assignment operator deleted
```

## Main Operators and Member Functions

### `AddEntry(const IndexType GlobalRowIndex, const IndexType GlobalColIndex)`
Adds a single directed edge (entry) from global row `GlobalRowIndex` to global column `GlobalColIndex`.
- If `GlobalRowIndex` is owned by the current rank (i.e., `GetRowNumbering().IsLocal(GlobalRowIndex)` is true), the entry `(LocalRowIndex, GlobalColIndex)` is added to the local graph part (`mLocalGraph`). `LocalRowIndex` is obtained via `GetRowNumbering().LocalId(GlobalRowIndex)`.
- If `GlobalRowIndex` is owned by another rank, the entry `(RemoteLocalRowIndex, GlobalColIndex)` is temporarily stored in a per-rank buffer (`mNonLocalGraphs[owner_rank]`). `RemoteLocalRowIndex` is obtained via `GetRowNumbering().RemoteLocalId(GlobalRowIndex, owner_rank)`. These entries are sent to their owner ranks during `Finalize()`.
A lock is used for thread-safe addition to non-local buffers.
```cpp
// Assuming dist_graph is initialized.
Kratos::IndexType global_row_idx = 5;  // Example global row index
Kratos::IndexType global_col_idx = 10; // Example global column index

dist_graph.AddEntry(global_row_idx, global_col_idx);
// If current rank owns global_row_idx 5, (5,10) is added to its local graph.
// Otherwise, it's stored to be sent to the owner of row 5 later.
```

### `AddEntries(const IndexType GlobalRowIndex, const std::vector<IndexType>& rGlobalColIndices)`
Adds multiple directed edges from a source global row `GlobalRowIndex` to all global column indices specified in the container `rGlobalColIndices`.
The handling of local vs. non-local `GlobalRowIndex` (including conversion to local/remote-local IDs and use of locks) is the same as for `AddEntry`. Column indices in `rGlobalColIndices` are global.
```cpp
Kratos::IndexType source_row = 7;
std::vector<Kratos::IndexType> target_cols = {12, 15, 20};

dist_graph.AddEntries(source_row, target_cols);
// Adds (7,12), (7,15), (7,20), handled based on ownership of row 7.
```
### `AddEntries(const std::vector<IndexType>& rGlobalRowIndices, const std::vector<IndexType>& rGlobalColIndices)`
For each global row index `I` in the container `rGlobalRowIndices`, adds directed edges from `I` to all global column indices in the container `rGlobalColIndices`.
This is equivalent to calling `AddEntries(I, rGlobalColIndices)` for each `I` in `rGlobalRowIndices`.
```cpp
std::vector<Kratos::IndexType> source_rows = {0, 4};
std::vector<Kratos::IndexType> target_cols_for_multiple = {1, 3};
dist_graph.AddEntries(source_rows, target_cols_for_multiple);
// Adds (0,1), (0,3) and (4,1), (4,3), handled based on ownership of rows 0 and 4.
```

### `AddEntries(const std::vector<std::vector<IndexType>>& rGlobalConnectivities)`
Adds entries based on a list of connectivities. Each inner vector (representing, e.g., an element's node global IDs) in `rGlobalConnectivities` defines a set of rows for which entries will be added.
For each global row index `I` within an inner vector `connectivities_I_J_K`, this function adds directed edges from `I` to all global column indices also present in `connectivities_I_J_K` (i.e., `J` and `K`, and `I` itself).
The handling of local vs. non-local rows (including conversion to local/remote-local IDs and use of locks) is the same as for `AddEntry`.
This is typically used to add all connections within an element.
```cpp
std::vector<std::vector<Kratos::IndexType>> connectivities;
connectivities.push_back({0,1,2}); // Element 1 connecting global nodes 0,1,2
connectivities.push_back({2,3,0}); // Element 2 connecting global nodes 2,3,0

dist_graph.AddEntries(connectivities);
// For {0,1,2}: attempts to add (0,0),(0,1),(0,2), (1,0),(1,1),(1,2), (2,0),(2,1),(2,2).
// Ownership of rows 0, 1, 2 determines where these are initially stored.
```
### `Finalize()`
Finalizes the graph construction. This is a **collective operation** that must be called on all MPI ranks participating in the `DataCommunicator` of the graph.
It performs crucial communication steps:
1. Exchanges graph parts: Each rank sends the non-local entries it has collected (from its `mNonLocalGraphs` buffers) to the MPI ranks that own those rows.
2. Merges received data: Each rank receives graph entries for rows it owns (sent by other ranks) and merges them into its `mLocalGraph`. This merge is done by converting the received `SingleVectorRepresentation` of a `SparseGraph` using `mLocalGraph.AddFromSingleVectorRepresentation()`.
The `mNonLocalGraphs` buffers are effectively send buffers and are cleared after their content is processed or sent, although the `Finalize` implementation does not explicitly show them being cleared (they are stateful until next `Clear()` or destruction).
The primary outcome of `Finalize` is that `mLocalGraph` on each rank contains all entries for its owned rows, whether they originated locally or were received from other ranks.
```cpp
dist_graph.Finalize(); // Must be called by all ranks.
```

### `GetLocalGraph() const`
Returns a constant reference to the local part of the graph stored on the current MPI rank. This graph contains entries `(LocalRowIndex, GlobalColIndex)` where `LocalRowIndex` is the local ID of a row owned by the current rank, and `GlobalColIndex` are global column indices.
The returned type is `const SparseContiguousRowGraph<TIndexType>&`.
```cpp
const Kratos::SparseContiguousRowGraph<Kratos::IndexType>& local_graph = dist_graph.GetLocalGraph();

// To iterate, using local row indices:
for(Kratos::IndexType local_row_idx = 0; local_row_idx < local_graph.Size(); ++local_row_idx) {
    // To get global row index:
    // Kratos::IndexType global_row_idx = dist_graph.GetRowNumbering().GlobalId(local_row_idx);
    const auto& neighbors = local_graph[local_row_idx]; // Neighbors are global column indices
    // Process local_row_idx (or global_row_idx) and its neighbors...
}
```

### `GetNonLocalGraph(IndexType Rank) const`
Returns a constant reference to a `SparseGraph<TIndexType>` that contains entries `(RemoteLocalRowIndex, GlobalColIndex)` which this rank has buffered to send to the specified `Rank`. `RemoteLocalRowIndex` is the local ID of the row on the target `Rank`.
This is primarily for internal inspection or debugging as these buffers are processed by `Finalize()`.
```cpp
Kratos::IndexType target_rank = 1; // Example rank
const Kratos::SparseGraph<Kratos::IndexType>& non_local_buffer_for_rank_1 = dist_graph.GetNonLocalGraph(target_rank);
// This graph contains entries this rank wants to add to rows owned by rank 1.
```

### `GetNonLocalGraphs() const`
Returns a constant reference to the `DenseVector` of `SparseGraph<TIndexType>` objects. Each element `[i]` in this vector is the buffer of non-local entries destined for rank `i`.
```cpp
const Kratos::DenseVector<Kratos::SparseGraph<TIndexType>>& all_non_local_buffers = dist_graph.GetNonLocalGraphs();
// all_non_local_buffers[some_rank] is the same as dist_graph.GetNonLocalGraph(some_rank).
```

### `GetRowNumbering() const`
Returns a constant reference to the internally managed `DistributedNumbering<TIndexType>` object. This object handles the distribution of rows and the mapping between local and global indices.
```cpp
const Kratos::DistributedNumbering<Kratos::IndexType>& numbering = dist_graph.GetRowNumbering();
// Kratos::IndexType owner_rank = numbering.OwnerRank(some_global_row_index);
// bool is_local = numbering.IsLocal(some_global_row_index);
```

### `Size()`
Returns the total number of global row indices managed by the `DistributedNumbering` across all ranks. This is the global dimension (number of rows) of the conceptual distributed matrix.
```cpp
Kratos::IndexType global_number_of_rows = dist_graph.Size();
// This value is determined by the DistributedNumbering setup.
```

### `LocalSize()`
Returns the number of row indices owned by the current MPI rank, as managed by the `DistributedNumbering`.
This is equivalent to `dist_graph.GetRowNumbering().LocalSize()`.
```cpp
Kratos::IndexType num_locally_owned_rows_with_entries = dist_graph.LocalSize();
```

### `Has(const IndexType GlobalRowIndex, const IndexType GlobalColIndex) const`
Checks if a specific edge `(GlobalRowIndex, GlobalColIndex)` exists in the **local part** of the graph.
The `GlobalRowIndex` is first converted to its corresponding `LocalRowIndex` using `GetRowNumbering().LocalId(GlobalRowIndex)`. The check is then performed as `mLocalGraph.Has(LocalRowIndex, GlobalColIndex)`.
This function will only yield meaningful results (and avoid errors) if `GlobalRowIndex` is actually owned by the current rank.
```cpp
Kratos::IndexType global_row_idx = 5;
Kratos::IndexType global_col_idx = 10;

bool entry_exists_locally = false;
if (dist_graph.GetRowNumbering().IsLocal(global_row_idx)) {
   entry_exists_locally = dist_graph.Has(global_row_idx, global_col_idx);
}

if (entry_exists_locally) {
    // Entry (5,10) exists and row 5 is owned by this rank.
}
```

### `Clear()`
Clears all data in the `DistributedSparseGraph` on the current rank. This includes clearing `mLocalGraph` (the local part of the graph) and `mNonLocalGraphs` (the send buffers for non-local entries).
This operation is local to the calling MPI rank. To clear the entire distributed graph, it must be called on all ranks.
```cpp
dist_graph.Clear();
// Local data structures are cleared.
// To reset the entire distributed graph, call Clear() on all ranks.
```

### `operator[](const IndexType& LocalRowIndex) const`
Provides access to the neighbors of a locally owned row, using its `LocalRowIndex`.
Returns a constant reference to the set of global column indices for that row.
The type of the returned object is `const typename LocalGraphType::GraphType::value_type&`, which for `SparseContiguousRowGraph` typically means `const std::vector<IndexType>&` or similar, representing the column entries for `LocalRowIndex`.
```cpp
// Assuming dist_graph is finalized and populated.
Kratos::IndexType local_row_to_access = 0; // Example: Access first local row
if (local_row_to_access < dist_graph.GetLocalGraph().Size()) {
    const auto& neighbors = dist_graph[local_row_to_access];
    std::cout << "Neighbors of local row " << local_row_to_access << ": ";
    for (Kratos::IndexType global_col_idx : neighbors) {
        std::cout << global_col_idx << " ";
    }
    std::cout << std::endl;
}
```

### `ComputeLocalMinMaxColumnIndex(IndexType& rMinJ, IndexType& rMaxJ) const`
Computes the minimum and maximum global column indices (`rMinJ`, `rMaxJ`) present in the entries of the local graph part (`mLocalGraph`).
```cpp
Kratos::IndexType local_min_col, local_max_col;
dist_graph.ComputeLocalMinMaxColumnIndex(local_min_col, local_max_col);
// local_min_col and local_max_col now hold the min/max global column indices on this rank.
```

### `ComputeMaxGlobalColumnIndex() const`
Computes the maximum global column index across all entries in the entire distributed graph. This involves finding the local maximum on each rank (using `ComputeLocalMinMaxColumnIndex`) and then performing a global `MaxAll` reduction.
Returns the global maximum column index.
```cpp
Kratos::IndexType max_global_j = dist_graph.ComputeMaxGlobalColumnIndex();
// max_global_j is the largest column index present anywhere in the distributed graph.
```

## Input and Output Methods

### `Info() const`
Returns a string containing basic information about the `DistributedSparseGraph`. As per the header, it returns the string "DistributedSparseGraph".
```cpp
std::cout << dist_graph.Info() << std::endl;
// Output: DistributedSparseGraph
```

### `PrintInfo(std::ostream& rOStream) const`
Prints basic information about the `DistributedSparseGraph` to the specified output stream `rOStream`. It prints the string "DistributedSparseGraph".
This method itself is not collective, it prints from the calling rank. To print only from rank 0:
```cpp
if (dist_graph.GetComm().Rank() == 0) {
    dist_graph.PrintInfo(std::cout); // Prints "DistributedSparseGraph"
    std::cout << std::endl;
}
```

### `PrintData(std::ostream& rOStream) const`
Prints the data of the `DistributedSparseGraph`.
**Note:** The default implementation of `PrintData` in `distributed_sparse_graph.h` is empty (`{}`), so it prints nothing.
To visualize the graph's contents, one would need to iterate through `GetLocalGraph()` on each rank and print manually, coordinating output if printing to a shared stream.
```cpp
Kratos::DataCommunicator& comm = dist_graph.GetRowNumbering().GetDataCommunicator();
for (int i = 0; i < comm.Size(); ++i) {
    comm.Barrier();
    if (comm.Rank() == i) {
        std::cout << "Data for rank " << i << ":" << std::endl;
        dist_graph.PrintData(std::cout); // Prints local graph part
        std::cout << std::endl;
    }
}
```

### Stream Operators
**`operator<<(std::ostream& rOStream, const DistributedSparseGraph<TIndexType>& rThis)`**
Overloads the `<<` operator for output streams. It typically calls `PrintInfo` (which prints on rank 0) and then `PrintData` (which prints local data from each rank, potentially leading to interleaved output if not managed).
```cpp
// If printing to a shared stream like std::cout from multiple ranks:
// On rank 0:
if (dist_graph.GetRowNumbering().GetDataCommunicator().Rank() == 0) {
   std::cout << "Graph Info (from rank 0):\n" << dist_graph; // Will call PrintInfo
}
// if (dist_graph.GetComm().Rank() == 0) { // To print only from rank 0
//    std::cout << dist_graph << std::endl;
// }
// Output (if rank 0 prints, and PrintData is empty):
// DistributedSparseGraph
//
```

## Serialization

The `DistributedSparseGraph` class provides `save` and `load` methods for serialization.
**Note:** As of the current header file (`kratos/containers/distributed_sparse_graph.h`), these methods are implemented as empty stubs:
`void save(Serializer& rSerializer) const {}`
`void load(Serializer& rSerializer) {}`
Therefore, direct serialization using these methods will not save or load any data for `DistributedSparseGraph`. If serialization is required, custom logic or a different mechanism would be needed, or these methods would need to be implemented fully.

### `save(Serializer& rSerializer) const` (Stub)
This method is intended to serialize the `DistributedSparseGraph` object. In the current implementation, it is an empty function.
```cpp
// dist_graph.save(rSerializer); // Calling this will do nothing.
```

### `load(Serializer& rSerializer)` (Stub)
This method is intended to deserialize the `DistributedSparseGraph` object. In the current implementation, it is an empty function.
```cpp
// dist_graph.load(rSerializer); // Calling this will do nothing.
```