---
title: SparseGraph
keywords: SparseGraph, Sparse Matrix, CSR
tags: [SparseGraph, Data Structures, Linear Algebra]
sidebar: kratos_for_developers
summary: Comprehensive documentation for the Kratos SparseGraph class, a sparse graph.
---

## SparseGraph Class

The `SparseGraph` class implements a sparse graph data structure. It is designed to store a matrix graph, aimed at the **fast construction** of other sparse matrix formats (particularly CSR).
It efficiently stores and manages graph connectivity where the number of edges is significantly smaller than the number of possible edges. This class is particularly useful in finite element analysis and other numerical simulations where sparse matrices are common.

**Important Note:** This class is **not thread-safe by design**. A graph should be computed in each thread and then merged if parallel computation is needed.

### Template Parameters

The `SparseGraph` class has one template parameter:

- `TIndexType`: Specifies the type of indices used for the graph nodes. This is typically an integer type, such as `std::size_t`.

### Lifecycle Methods

#### Constructors

**Default Constructor**
Creates an empty `SparseGraph`.
```cpp
Kratos::SparseGraph<std::size_t> graph;
```

**Constructor with Size Hint `N`**
`SparseGraph(IndexType N)`: This constructor initializes an empty graph. The parameter `N` is a hint and is not used to pre-allocate memory for the graph structure itself (which is a `std::map` and grows dynamically). The `Size()` of the graph is determined by the maximum row index to which entries are added, not by `N` directly.
```cpp
Kratos::SparseGraph<std::size_t> graph(10); // N=10 is passed, graph is initially empty.
                                            // Its Size() will be 0 until entries are added.
```
**Constructor with DataCommunicator**
`SparseGraph(DataCommunicator& rComm)`: Creates a `SparseGraph` associated with a specific Kratos `DataCommunicator`.
An error is raised if `rComm.IsDistributed()` is true, as `SparseGraph` is designed for serial execution.
```cpp
Kratos::DataCommunicator& r_comm = Kratos::ParallelEnvironment::GetDefaultDataCommunicator();
Kratos::SparseGraph<std::size_t> graph(r_comm);
```

**Constructor from Single Vector Representation**
`SparseGraph(const std::vector<IndexType>& rSingleVectorRepresentation)`: Creates a `SparseGraph` by populating it from the data provided in `rSingleVectorRepresentation`. This vector follows a specific format described in `ExportSingleVectorRepresentation` and `AddFromSingleVectorRepresentation`.
```cpp
// Example: represents a graph with 3 rows (max index 2): 0->{1,2}, 2->{0}
std::vector<std::size_t> single_vector_rep = {3, 0, 2, 1, 2, 2, 1, 0};
Kratos::SparseGraph<std::size_t> graph_from_vec(single_vector_rep);
```
#### Destructor
The destructor handles the deallocation of resources managed by the `SparseGraph`.
```cpp
// graph object goes out of scope, its destructor is automatically called.
```

#### Copy Constructor
`SparseGraph(const SparseGraph& rOther)`: Creates a new `SparseGraph` as a deep copy of `rOther`.
```cpp
Kratos::SparseGraph<std::size_t> original_graph;
original_graph.AddEntry(0,1);
// ... populate original_graph ...
Kratos::SparseGraph<std::size_t> copied_graph(original_graph);
// copied_graph is now an independent copy of original_graph.
```
#### Assignment Operator
`SparseGraph& operator=(SparseGraph const& rOther) = delete;`
The assignment operator is explicitly deleted. `SparseGraph` objects cannot be assigned to one another after construction.
```cpp
Kratos::SparseGraph<std::size_t> graph1, graph2;
graph1 = graph2; // This line would cause a compilation error.
```

### Main Operators and Member Functions

#### `AddEntry(TIndexType Index1, TIndexType Index2)`
Adds a single directed edge (entry) from `Index1` to `Index2`. If the entry already exists, it is not added again.
```cpp
Kratos::SparseGraph<std::size_t> graph(5);
graph.AddEntry(0, 1); // Adds an edge from node 0 to node 1
graph.AddEntry(0, 2);
```
#### `AddEntries(TIndexType Index, const std::vector<TIndexType>& rNewEntries)`
Adds multiple directed edges from a source node `Index` to all nodes specified in `rNewEntries`.
```cpp
Kratos::SparseGraph<std::size_t> graph(5);
std::vector<std::size_t> targets = {1, 2, 3};
graph.AddEntries(0, targets); // Adds edges from node 0 to nodes 1, 2, and 3
```

#### `AddEntries`
There are several overloads to add multiple entries:

**`AddEntries(const IndexType RowIndex, const TContainerType& rColIndices)`**
Adds multiple directed edges from a source node `RowIndex` to all nodes specified in the container `rColIndices`.
```cpp
Kratos::SparseGraph<std::size_t> graph;
std::vector<std::size_t> targets = {1, 2, 3};
graph.AddEntries(0, targets); // Adds edges from node 0 to nodes 1, 2, and 3
```
**`AddEntries(const IndexType RowIndex, const TIteratorType& rColBegin, const TIteratorType& rColEnd)`**
Adds multiple directed edges from `RowIndex` to column indices specified by the iterator range `[rColBegin, rColEnd)`.
```cpp
Kratos::SparseGraph<std::size_t> graph;
std::vector<std::size_t> targets = {4, 5};
graph.AddEntries(0, targets.begin(), targets.end()); // Adds edges from node 0 to nodes 4 and 5
```

**`AddEntries(const TContainerType& rIndices)`**
For each index `I` in `rIndices`, adds directed edges from `I` to all other indices in `rIndices`.
Effectively, for each node `I` in `rIndices`, entries `(I, J)` are added for every node `J` in `rIndices`.
```cpp
Kratos::SparseGraph<std::size_t> graph;
std::vector<std::size_t> nodes_to_connect = {0, 1, 2};
graph.AddEntries(nodes_to_connect);
// This results in:
// Row 0 gets neighbors {0,1,2}
// Row 1 gets neighbors {0,1,2}
// Row 2 gets neighbors {0,1,2}
// (Order within each std::unordered_set is not guaranteed).
```
**`AddEntries(const TContainerType& rRowIndices, const TContainerType& rColIndices)`**
For each index `I` in `rRowIndices`, adds directed edges from `I` to all indices in `rColIndices`.
```cpp
Kratos::SparseGraph<std::size_t> graph;
std::vector<std::size_t> sources = {0, 1};
std::vector<std::size_t> destinations = {2, 3};
graph.AddEntries(sources, destinations);
// Resulting graph: 0->{2,3}, 1->{2,3}
```

**`AddEntries(SparseGraph& rOtherGraph)`**
Merges all entries from `rOtherGraph` into the current graph.
```cpp
Kratos::SparseGraph<std::size_t> graph1, graph2;
graph1.AddEntry(0,1);
graph2.AddEntry(0,2);
graph2.AddEntry(1,2);
graph1.AddEntries(graph2);
// graph1 now contains: 0->{1,2}, 1->{2}
```

#### `ExportCSRArrays`
Exports the graph connectivity in Compressed Sparse Row (CSR) format. The column indices for each row are sorted internally by this function.

**`IndexType ExportCSRArrays(TVectorType& rRowIndices, TVectorType& rColIndices) const`**
Exports the graph to provided vector containers.
- `rRowIndices`: Will be resized and filled with row pointers. Its size will be `nrows + 1`, where `nrows` is the number of rows in the graph (i.e., `this->Size()`).
- `rColIndices`: Will be resized and filled with column indices. Its size will be `nnz` (the total number of non-zero entries).
Returns `nrows` (the number of rows, equal to `this->Size()`).
```cpp
Kratos::SparseGraph<std::size_t> graph;
graph.AddEntry(0,1);
graph.AddEntry(0,2); // graph: 0 -> {1,2}
graph.AddEntry(1,2); // graph: 1 -> {2}
// Max row index is 1. So, Size() will be 2.

std::vector<std::size_t> row_indices_vec;
std::vector<std::size_t> col_indices_vec;
std::size_t num_rows = graph.ExportCSRArrays(row_indices_vec, col_indices_vec);

// num_rows will be 2.
// row_indices_vec will be [0, 2, 3]
// (Row 0 has 2 entries, Row 1 has 1 entry)
// col_indices_vec will be [1, 2, 2]
// (Entries for row 0 are {1,2}, for row 1 is {2}. Columns are sorted per row.)
```

**`IndexType ExportCSRArrays(IndexType*& pRowIndicesData, IndexType& rRowDataSize, IndexType*& pColIndicesData, IndexType& rColDataSize) const`**
This version allocates memory for `pRowIndicesData` (size `nrows+1`) and `pColIndicesData` (size `nnz`) using `new[]`.
The caller is responsible for `delete[]`ing these arrays.
- `rRowDataSize`: Output, set to `nrows + 1`.
- `rColDataSize`: Output, set to `nnz`.
Returns `nrows`.
```cpp
Kratos::SparseGraph<std::size_t> graph;
graph.AddEntry(0,1); // Example graph
graph.AddEntry(2,0); // Max row index 2, so nrows = 3

std::size_t* p_row_indices = nullptr;
std::size_t row_size = 0;
std::size_t* p_col_indices = nullptr;
std::size_t col_size = 0;

std::size_t num_rows_exported = graph.ExportCSRArrays(p_row_indices, row_size, p_col_indices, col_size);
// num_rows_exported = 3
// row_size = 4 (nrows+1)
// p_row_indices might be [0, 1, 1, 2] (0->{1}, 1->{}, 2->{0})
// col_size = 2 (total entries)
// p_col_indices might be [1, 0] (cols for row 0, then cols for row 2)

delete[] p_row_indices;
delete[] p_col_indices;
```

**`IndexType ExportCSRArrays(Kratos::span<IndexType>& rRowIndices, Kratos::span<IndexType>& rColIndices) const = delete;`**
This overload using `Kratos::span` is explicitly deleted and cannot be used.

#### `ExportSingleVectorRepresentation() const`
Exports the graph connectivity as a single `std::vector<IndexType>`. This vector is returned by the function.
The format of the returned vector is:
`[NROWS, ROW_0_INDEX, R0_NUM_NEIGHBORS, R0_neighbor1, R0_neighbor2, ..., ROW_I_INDEX, RI_NUM_NEIGHBORS, RI_neighbor1, ...]`
- The first element `NROWS` is the total number of rows in the graph, as given by `this->Size()`.
- This is followed by blocks of data for each non-empty row in the graph. Each block consists of:
    - The row index.
    - The number of neighbors (columns) for that row.
    - The column indices of the neighbors.
```cpp
Kratos::SparseGraph<std::size_t> graph;
graph.AddEntry(0,1);
graph.AddEntry(0,2); // graph: 0 -> {1,2} (max index 0, so Size is 1, but if we add 2,0 then Size is 3)
graph.AddEntry(2,0); // graph: 2 -> {0} (max index 2, so Size is 3)

std::vector<std::size_t> single_vec_rep = graph.ExportSingleVectorRepresentation();
// If graph is 0->{1,2}, 2->{0}:
// Size() is 3 (max row index 2 + 1).
// Row 0: index 0, num_neighbors 2, neighbors {1,2}.
// Row 1: is empty, so it's not explicitly listed in the data part.
// Row 2: index 2, num_neighbors 1, neighbor {0}.
// single_vec_rep would be: [3, 0, 2, 1, 2, 2, 1, 0]
// (Order of neighbors like {1,2} for row 0 is not guaranteed as it comes from an unordered_set).
```

#### `AddFromSingleVectorRepresentation(const std::vector<IndexType>& rSingleVectorRepresentation)`
Clears the current graph and populates it from the data provided in `rSingleVectorRepresentation`, which must follow the format described in `ExportSingleVectorRepresentation`.
```cpp
Kratos::SparseGraph<std::size_t> graph;
std::vector<std::size_t> rep_data = {3, 0, 2, 1, 2, 2, 1, 0}; // Represents 0->{1,2}, 2->{0}
graph.AddFromSingleVectorRepresentation(rep_data);
// graph now contains entries (0,1), (0,2), and (2,0).
```
#### `GetGraph() const`
Returns a constant reference to the underlying graph data structure.
The graph is stored as `std::map<IndexType, std::unordered_set<IndexType> >`.
```cpp
Kratos::SparseGraph<std::size_t> graph;
graph.AddEntry(0,1);
graph.AddEntry(0,2);
const std::map<std::size_t, std::unordered_set<std::size_t>>& internal_graph_data = graph.GetGraph();
for(const auto& row_pair : internal_graph_data) {
    std::cout << "Row " << row_pair.first << " has neighbors: ";
    for(std::size_t col_idx : row_pair.second) {
        std::cout << col_idx << " ";
    }
    std::cout << std::endl;
}
// Output for graph 0->{1,2}:
// Row 0 has neighbors: 1 2  (or 2 1, as std::unordered_set does not guarantee order)
```

#### `Size() const`
Returns the "size" of the graph. This is calculated as the largest row index that has entries, plus one. If the graph is empty (no entries), it returns 0.
```cpp
Kratos::SparseGraph<std::size_t> graph;
std::cout << "Initial size: " << graph.Size() << std::endl; // Output: Initial size: 0
graph.AddEntry(2, 0); // Max row index is 2
std::cout << "Size after adding (2,0): " << graph.Size() << std::endl; // Output: Size after adding (2,0): 3
graph.AddEntry(0, 1); // Max row index is still 2
std::cout << "Size after adding (0,1): " << graph.Size() << std::endl; // Output: Size after adding (0,1): 3
```

#### `IsEmpty() const`
Returns `true` if the graph contains no entries (i.e., `mGraph` is empty), `false` otherwise.
```cpp
Kratos::SparseGraph<std::size_t> graph;
KRATOS_CHECK(graph.IsEmpty()); // True
graph.AddEntry(0,0);
KRATOS_CHECK_IS_FALSE(graph.IsEmpty()); // False
```
#### `Has(TIndexType Index1, TIndexType Index2)`
Checks if a directed edge exists from `Index1` to `Index2`.
Returns `true` if the edge exists, `false` otherwise.
```cpp
Kratos::SparseGraph<std::size_t> graph(3);
graph.AddEntry(0, 1);
bool exists = graph.Has(0, 1); // true
bool not_exists = graph.Has(1, 0); // false (unless explicitly added)
```

#### `operator[](const IndexType& Key) const`
Returns a constant reference to the `std::unordered_set<IndexType>` of neighbors for the row `Key`.
**Important**: This operator assumes that `Key` exists as a row in the graph (i.e., `mGraph.find(Key)` would not be `mGraph.end()`). Accessing a non-existent key via this `operator[]` can lead to undefined behavior because it would involve dereferencing an end iterator of the underlying map.
It is safer to first check if the row `Key` exists using `graph.GetGraph().count(Key)` or by iterating.
```cpp
Kratos::SparseGraph<std::size_t> graph;
graph.AddEntry(0, 1);
graph.AddEntry(0, 2);

// Safe access:
if (graph.GetGraph().count(0)) { // Check if row 0 exists
    const std::unordered_set<std::size_t>& neighbors_of_0 = graph[0];
    std::cout << "Neighbors of 0: ";
    for (std::size_t neighbor : neighbors_of_0) {
        std::cout << neighbor << " "; // Output: 1 2 (or 2 1, order not guaranteed)
    }
    std::cout << std::endl;
}

// Accessing graph[1] here would be unsafe if row 1 has not been added.
const auto& neighbors_of_1 = graph[1]; // Potential Undefined Behavior if 1 is not a key in mGraph
```

#### `Clear()`
Removes all nodes and edges from the graph, making it empty.
```cpp
Kratos::SparseGraph<std::size_t> graph(5);
graph.AddEntry(0,1);
graph.Clear();
std::cout << "Graph size after clear: " << graph.Size() << std::endl; // Output: Graph size after clear: 0
```
#### Iterators
The `SparseGraph` provides iterators to traverse the graph data.
- `begin()`: Returns an iterator to the beginning of the graph's adjacency list.
- `end()`: Returns an iterator to the end of the graph's adjacency list.
- `cbegin()`: Returns a const iterator to the beginning.
- `cend()`: Returns a const iterator to the end.
```cpp
Kratos::SparseGraph<std::size_t> graph(3);
graph.AddEntry(0, 1);
graph.AddEntry(1, 2);
graph.AddEntry(0, 2);
The `SparseGraph` provides a `const_iterator_adaptor` to iterate over its rows (the entries in the underlying `std::map`).
- `begin()`: Returns a `const_iterator_adaptor` to the first row of the graph.
- `end()`: Returns a `const_iterator_adaptor` to the end of the graph entries.
The iterator dereferences to `const typename GraphType::mapped_type&`, which is `const std::unordered_set<IndexType>&` (the set of column indices for the current row).
The iterator also has a `GetRowIndex()` method to get the current row index.
```cpp
Kratos::SparseGraph<std::size_t> graph;
graph.AddEntry(0, 1);
graph.AddEntry(1, 2);
graph.AddEntry(0, 2); // Row 0: {1,2}, Row 1: {2}
std::cout << "Graph iteration using iterators:" << std::endl;
for (auto it_row = graph.begin(); it_row != graph.end(); ++it_row) {
    std::cout << "Row " << it_row.GetRowIndex() << " has neighbors: ";
    const std::unordered_set<std::size_t>& neighbors = *it_row; // Dereference iterator
    for (const auto& neighbor_idx : neighbors) {
        std::cout << neighbor_idx << " ";
    }
    std::cout << std::endl;
}
// Example Output:
// Row 0 has neighbors: 1 2  (or 2 1)
// Row 1 has neighbors: 2 
```
#### `Finalize()`
This method is provided but currently has an empty implementation in `sparse_graph.h`. It might be used in derived classes or future versions for graph optimization or finalization steps.
```cpp
Kratos::SparseGraph<std::size_t> graph;
graph.AddEntry(0,1);
graph.Finalize(); // Currently does nothing for SparseGraph itself.
```

#### `GetComm() const` and `pGetComm() const`
Return the `DataCommunicator` associated with the graph. `SparseGraph` is intended for serial use, so this will typically be the serial communicator.
- `GetComm()`: Returns a const reference to the communicator.
- `pGetComm()`: Returns a const pointer to the communicator.
```cpp
Kratos::SparseGraph<std::size_t> graph;
const Kratos::DataCommunicator& comm_ref = graph.GetComm();
// const Kratos::DataCommunicator* comm_ptr = graph.pGetComm();
KRATOS_CHECK_IS_FALSE(comm_ref.IsDistributed()); // Should be true for default serial communicator
```

### Input and Output Methods

#### `Info() const`
Returns a string containing basic information about the `SparseGraph`. According to `sparse_graph.h`, this returns the string "SparseGraph".
```cpp
Kratos::SparseGraph<std::size_t> graph;
std::cout << graph.Info() << std::endl;
// Output: SparseGraph
```

#### `PrintInfo(std::ostream& rOStream) const`
Prints basic information about the `SparseGraph` to the specified output stream `rOStream`. It prints the string "SparseGraph".
```cpp
Kratos::SparseGraph<std::size_t> graph;
graph.PrintInfo(std::cout); // Prints "SparseGraph" to cout
std::cout << std::endl;
```

#### `PrintData(std::ostream& rOStream) const`
Prints the data of the `SparseGraph` to the specified output stream `rOStream`.
**Note:** The default implementation of `PrintData` in `sparse_graph.h` is empty, so it prints nothing.
```cpp
Kratos::SparseGraph<std::size_t> graph;
graph.AddEntry(0,1);
graph.PrintData(std::cout); // This line will print nothing.
std::cout << std::endl;
```
To visualize the graph's contents, iterate using `GetGraph()` or the graph's iterators and print manually.

#### Stream Operators

**`operator<<(std::ostream& rOStream, const SparseGraph<TIndexType>& rThis)`**
Overloads the `<<` operator for output streams. It calls `rThis.PrintInfo(rOStream)`, then prints a newline, and then calls `rThis.PrintData(rOStream)`.
```cpp
Kratos::SparseGraph<std::size_t> graph;
graph.AddEntry(0,1);
std::cout << graph;
// Example output (based on sparse_graph.h where PrintData is empty):
// SparseGraph
// 
// (The output is "SparseGraph" followed by a newline, then whatever PrintData outputs (nothing), 
//  then the stream operator itself doesn't add another newline automatically after PrintData)
```

### Serialization

#### `save(Serializer& rSerializer)`
Serializes the `SparseGraph` object. This is typically used for saving the state of the graph to a file or sending it over a network.
```cpp
// Assuming rSerializer is a Kratos::Serializer object
Kratos::SparseGraph<std::size_t> graph_to_save;
// ... populate graph_to_save ...
rSerializer.save("MyGraph", graph_to_save);
```
**Note:** Actual usage depends on the Kratos serialization framework.

#### `load(Serializer& rSerializer)`
Deserializes the `SparseGraph` object, restoring its state from a serialized representation.
```cpp
// Assuming rSerializer is a Kratos::Serializer object
Kratos::SparseGraph<std::size_t> graph_to_load;
rSerializer.load("MyGraph", graph_to_load);
```