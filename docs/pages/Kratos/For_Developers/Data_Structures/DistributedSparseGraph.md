---
title: DistributedSparseGraph
keywords: DistributedSparseGraph, Sparse Matrix, CSR, Distributed Computing, MPI, Parallel Linear Algebra
tags: [DistributedSparseGraph, Data Structures, Linear Algebra, Distributed, MPI, HPC]
sidebar: kratos_for_developers
summary: Comprehensive documentation for the Kratos DistributedSparseGraph class, a distributed sparse graph designed for high-performance parallel computations using MPI.
---

## DistributedSparseGraph Class

The `DistributedSparseGraph` class implements a distributed sparse graph data structure designed for parallel computations using **MPI**. It allows different **MPI** ranks to manage portions of a global graph, facilitating distributed assembly and storage of graph connectivity. This class is crucial for large-scale simulations where a single machine's memory is insufficient to hold the entire graph.

The graph is built upon a `DataCommunicator` which defines the **MPI** context. It internally creates and manages a `DistributedNumbering` object to handle the mapping between global and local row indices across ranks. Each rank stores its local portion of the graph (for rows it owns) and uses temporary storage for entries that belong to non-local rows, which are communicated during `Finalize()`.

**Important Note on Thread Safety:** While `DistributedSparseGraph` coordinates data across **MPI** ranks, operations that modify internal graph structures (especially those involving non-local entries before `Finalize()`) use locks for thread safety. However, the underlying local graph (`SparseContiguousRowGraph`) and non-local graph buffers (`SparseGraph`) might have their own considerations if populated directly in a multi-threaded fashion outside of the provided `AddEntry`/`AddEntries` methods. The header of the code itself notes: "it is BY DESIGN NOT threadsafe! (a graph should be computed in each thread and then merged)" - this generally applies to how graph data might be assembled before being passed to this class or if manipulating its internal structures directly.

For serial version see [SparseGraph](SparseGraph.md).

### Template Parameters

The `DistributedSparseGraph` class has one template parameter:

- `TIndexType`: Specifies the type of indices used for the graph nodes (e.g., `std::size_t`). This type should be consistent with the `DistributedNumbering` used.

### Lifecycle Methods

#### Constructor

**`DistributedSparseGraph(DistributedNumbering<TIndexType>& rRowNumbering)`**
Constructs a `DistributedSparseGraph`.
- `rRowNumbering`: A reference to a `DistributedNumbering` object. This object defines how global row indices are distributed among **MPI** ranks and provides methods for converting between global and local indices. The `DataCommunicator` is obtained from this numbering.

```cpp
// Assuming r_comm is an MPI_Comm (e.g., MPI_COMM_WORLD)
Kratos::DataCommunicator& kratos_comm = Kratos::ParallelEnvironment::GetDataCommunicator("World");

// Example: Create a DistributedNumbering for rows
// Suppose total 100 rows, rank 0 owns 0-49, rank 1 owns 50-99 (simplified)
std::vector<IndexType> partition;
// ... fill partition based on rank ...
// For rank 0: partition = {50} (meaning it owns 50 rows, indices 0 to 49 if ghost layers are 0)
// For rank 1: partition = {50} (meaning it owns 50 rows, indices 50 to 99 if ghost layers are 0)
// Actual partitioning depends on how DistributedNumbering is constructed.
// Let's assume a simple case for rRowNumbering construction:
Kratos::DistributedNumbering<IndexType> row_numbering(kratos_comm, partition);

// Construct the DistributedSparseGraph
Kratos::DistributedSparseGraph<IndexType> dist_graph(row_numbering);
```
#### Destructor
The destructor handles the deallocation of resources managed by the `DistributedSparseGraph`.
```cpp
// dist_graph goes out of scope, its destructor is automatically called.
```
**Note:** `DistributedSparseGraph` does not explicitly define a copy constructor or assignment operator, implying default behavior or that they are not intended for typical use.

### Main Operators and Member Functions

#### `AddEntry(const IndexType GlobalRowIndex, const IndexType GlobalColIndex)`
Adds a single directed edge (entry) from global row `GlobalRowIndex` to global column `GlobalColIndex`.
If `GlobalRowIndex` is owned by the current rank, the entry `(GlobalRowIndex, GlobalColIndex)` is added to the local part of the graph.
If `GlobalRowIndex` is owned by another rank, this information is stored temporarily and communicated during `Finalize()`.
```cpp
// Assuming dist_graph is initialized and rRowNumbering is available
IndexType my_rank = kratos_comm.Rank();
IndexType global_row_idx = /* some global row index */;
IndexType global_col_idx = /* some global column index */;

if (row_numbering.IsLocal(global_row_idx)) { // Example: rank 0 owns global rows 0-4
    dist_graph.AddEntry(global_row_idx, global_col_idx + my_rank); // Example entry
}
// Example from test_distributed_sparse_graph.cpp:
// rank 0 adds entry (0,0) and (0,1)
// rank 1 adds entry (2,2) and (2,3)
dist_graph.AddEntry(0,0); // If rank 0
dist_graph.AddEntry(2,2); // If rank 1
```

#### `AddEntries(const IndexType GlobalRowIndex, const std::vector<IndexType>& rGlobalColIndices)`
Adds multiple directed edges from a source global row `GlobalRowIndex` to all global column indices specified in `rGlobalColIndices`.
Similar to `AddEntry`, handles local and non-local `GlobalRowIndex` appropriately.
```cpp
IndexType global_row_idx = /* ... */;
std::vector<IndexType> global_target_cols = { /* ... */ };
dist_graph.AddEntries(global_row_idx, global_target_cols);
Example from test_distributed_sparse_graph.cpp:
std::vector<int> indices_to_add = {0,1}; // connectivities for a specific element
dist_graph.AddEntries(indices_to_add[0], indices_to_add); // Adds (0,0) and (0,1) if indices_to_add[0] is 0
```

#### `AddEntries(const std::vector<IndexType>& rGlobalRowIndices, const std::vector<IndexType>& rGlobalColIndices)`
For each global row index `GR_i` in `rGlobalRowIndices`, adds directed edges from `GR_i` to all global column indices in `rGlobalColIndices`.
```cpp
std::vector<IndexType> global_source_rows = { /* ... */ };
std::vector<IndexType> global_target_cols = { /* ... */ };
dist_graph.AddEntries(global_source_rows, global_target_cols);
```

#### `AddEntries(const std::vector<std::vector<IndexType>>& rGlobalConnectivities)`
Adds entries based on a list of connectivities. Each inner vector in `rGlobalConnectivities` represents a clique of interconnected global indices (e.g., nodes of a finite element). For each index `I` in an inner vector, entries `(I,J)` are added for all other indices `J` in the same inner vector.
```cpp
std::vector<std::vector<IndexType>> connectivities;
connectivities.push_back({0,1,2}); // Element 1 connecting nodes 0,1,2
connectivities.push_back({2,3,0}); // Element 2 connecting nodes 2,3,0
dist_graph.AddEntries(connectivities);
// For {0,1,2}: adds (0,0),(0,1),(0,2), (1,0),(1,1),(1,2), (2,0),(2,1),(2,2)
// (behavior regarding local/non-local rows applies)
```

#### `Finalize()`
Finalizes the graph construction. This is a collective operation that must be called on all **MPI** ranks.
It performs the following key steps:
1. Exchanges information about non-local entries: If a rank tried to add an entry for a row it does not own, this information is sent to the owner rank.
2. Populates the `mNonLocalGraph` on each rank: This graph stores, for rows owned by the current rank, which other ranks also have these rows as non-local (ghost) rows and need to know about their connectivity.
```cpp
dist_graph.Finalize(); // Collective call
```

#### `GetLocalGraph() const`
Returns a constant reference to the local part of the graph. This graph contains entries `(GlobalRowIndex, GlobalColIndex)` where `GlobalRowIndex` is owned by the current **MPI** rank.
The returned type is `const SparseGraph<TIndexType>&`.
```cpp
const Kratos::SparseGraph<IndexType>& local_graph = dist_graph.GetLocalGraph();
// Iterate over local_graph to see locally owned rows and their connections
for(auto it_row = local_graph.begin(); it_row != local_graph.end(); ++it_row) {
    IndexType global_row_idx = it_row.GetRowIndex();
    const auto& neighbors = *it_row;
    // ...
}
```

#### `GetNonLocalGraph() const`
Returns a constant reference to the non-local graph structure.
This is a `std::map<IndexType, std::vector<int>>`. The map key is a `GlobalRowIndex` (owned by the current rank). The associated `std::vector<int>` contains the ranks of the **MPI** processes that have this `GlobalRowIndex` as a ghost/non-local row and therefore need to be informed about its connectivity.
This is primarily for internal use or advanced scenarios.
```cpp
const std::map<IndexType, std::vector<int>>& non_local_graph_info = dist_graph.GetNonLocalGraph();
// For a locally owned GlobalRowIndex, non_local_graph_info[GlobalRowIndex] lists ranks that ghost it.
```

#### `GetRowNumbering() const`
Returns a constant reference to the `DistributedNumbering<TIndexType>` object used by the graph.
```cpp
const Kratos::DistributedNumbering<IndexType>& numbering = dist_graph.GetRowNumbering();
IndexType owner_rank = numbering.OwnerRank(global_idx);
```

#### `Size()`
Returns the total number of unique global row indices that have at least one entry in the graph, across all **MPI** ranks. This implies the "height" or number of active rows in the conceptual global matrix graph.
```cpp
IndexType total_active_rows = dist_graph.Size();
```

#### `LocalSize()`
Returns the number of unique global row indices owned by the current **MPI** rank that have at least one entry in the graph.
```cpp
IndexType local_active_rows = dist_graph.LocalSize();
// This is equivalent to dist_graph.GetLocalGraph().Size() if GetLocalGraph().Size() counts non-empty rows.
// More precisely, it's the number of keys in the local graph part.
```

#### `Has(const IndexType GlobalRowIndex, const IndexType GlobalColIndex) const`
Checks if a specific edge `(GlobalRowIndex, GlobalColIndex)` exists in the graph.
This operation is local if `GlobalRowIndex` is owned by the current rank. If `GlobalRowIndex` is owned by another rank, it currently might only check the sender's knowledge before `Finalize`, or it might not be well-defined for remote rows post-`Finalize` without further communication (the header suggests it checks `mpLocalGraph`). It's safest to assume this is primarily for locally owned rows.
```cpp
IndexType global_row_idx = /* ... */;
IndexType global_col_idx = /* ... */;
if (dist_graph.GetRowNumbering().IsLocal(global_row_idx)) {
   bool exists = dist_graph.Has(global_row_idx, global_col_idx);
   // ...
}
```

#### `Clear()`
Clears all data in the `DistributedSparseGraph` on all ranks. This includes local entries, non-local communication buffers, and resets the finalized state.
This is a collective operation.
```cpp
dist_graph.Clear(); // Collective call
// dist_graph is now empty and ready to be rebuilt.
```

### Input and Output Methods

#### `Info() const`
Returns a string containing basic information about the `DistributedSparseGraph`.
```cpp
std::cout << dist_graph.Info() << std::endl;
// Example output might include information about local/global size.
```

#### `PrintInfo(std::ostream& rOStream) const`
Prints basic information about the `DistributedSparseGraph` to the specified output stream `rOStream`.
This is a collective operation and typically prints information only on the root rank (rank 0).
```cpp
// Only rank 0 prints:
if (dist_graph.GetRowNumbering().GetDataCommunicator().Rank() == 0) {
    dist_graph.PrintInfo(std::cout);
    std::cout << std::endl;
}
```

#### `PrintData(std::ostream& rOStream) const`
Prints the data of the `DistributedSparseGraph` to the specified output stream `rOStream`.
This is a collective operation. Each rank prints its local data.
```cpp
// Each rank prints its local portion:
std::cout << "Data on Rank " << dist_graph.GetRowNumbering().GetDataCommunicator().Rank() << ":" << std::endl;
dist_graph.PrintData(std::cout);
std::cout << std::endl;

For controlled printing (e.g. rank by rank):
for (int i = 0; i < kratos_comm.Size(); ++i) {
    kratos_comm.Barrier();
    if (kratos_comm.Rank() == i) {
        std::cout << "Data for rank " << i << ":" << std::endl;
        dist_graph.PrintData(std::cout);
        std::cout << std::endl;
    }
    kratos_comm.Barrier();
}
```

#### Stream Operators
**`operator<<(std::ostream& rOStream, const DistributedSparseGraph<TIndexType>& rThis)`**
Overloads the `<<` operator for output streams. It typically calls `PrintInfo` and `PrintData`.
The behavior regarding which ranks print can be similar to `PrintInfo` or `PrintData`.
```cpp
if (dist_graph.GetRowNumbering().GetDataCommunicator().Rank() == 0) {
   std::cout << dist_graph << std::endl;
}
```

### Serialization

The `DistributedSparseGraph` class provides `save` and `load` methods for serialization, allowing its state to be saved and restored. This is essential for checkpointing and restarting large simulations.

#### `save(Serializer& rSerializer) const`
Serializes the `DistributedSparseGraph` object. This is a collective operation. Each rank serializes its local data.
```cpp
// Assuming rSerializer is a Kratos::Serializer object configured for distributed saving.
dist_graph.save(rSerializer);
```

#### `load(Serializer& rSerializer)`
Deserializes the `DistributedSparseGraph` object. This is a collective operation. Each rank deserializes its local data. The `DistributedNumbering` must be set up compatibly before loading.
```cpp
Kratos::DistributedNumbering<IndexType> row_numbering_for_load(...);
Kratos::DistributedSparseGraph<IndexType> graph_to_load(row_numbering_for_load);
graph_to_load.load(rSerializer);
```