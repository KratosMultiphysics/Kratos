---
title: DistributedSystemVector
keywords: DistributedSystemVector, Distributed Computing, MPI, Parallel Linear Algebra
tags: [DistributedSystemVector, Data Structures, Linear Algebra, Distributed, MPI, HPC]
sidebar: kratos_for_developers
summary: Comprehensive documentation for the Kratos SystemVector class, a distributed system of equations vector class designed for high-performance parallel computations using MPI.
---

## DistributedSystemVector Class

The `DistributedSystemVector` is a data structure designed for handling distributed vectors in parallel computations, primarily utilizing MPI (Message Passing Interface). It is a cornerstone in large-scale Finite Element Method (FEM) simulations where system vectors (such as solution vectors or right-hand side vectors) are partitioned across multiple processes.

Each MPI process owns a local portion of the vector. The `DistributedSystemVector` manages this local data and provides mechanisms for interacting with non-local data (data owned by other processes) when necessary, particularly during the assembly phase or global reduction operations (like norm or dot product).

### Constructors

`DistributedSystemVector` offers constructors tailored for parallel environments:

*   **From `DistributedSparseGraph`:**
    ```cpp
    DistributedSystemVector(const DistributedSparseGraph<IndexType>& rGraph);
    ```
    This is a common constructor used when the vector's structure is derived from a `DistributedSparseGraph`. It initializes the vector's local data size based on the graph's local row count and pre-allocates space in `mNonLocalData` for entries corresponding to non-local rows that this process will need to contribute to during assembly. The `DataCommunicator` and `DistributedNumbering` are also inherited from the graph.

*   **From `DistributedNumbering`:**
    ```cpp
    DistributedSystemVector(const DistributedNumbering<IndexType>& rNumbering);
    ```
    Constructs a `DistributedSystemVector` using an existing `DistributedNumbering` object. The local data vector is sized according to `rNumbering.LocalSize()`.
    **Important:** If this constructor is used, the vector cannot be directly used for assembly involving non-local contributions unless `AddEntry` or `AddEntries` is called to populate the `mNonLocalData` map *before* `BeginAssemble()` is called. This is because `mNonLocalData` (which determines which non-local entries are communicated) is not pre-filled as it is with the `DistributedSparseGraph` constructor.

*   **Copy Constructor:**
    ```cpp
    explicit DistributedSystemVector(DistributedSystemVector const& rOther);
    ```
    Creates a deep copy of another `DistributedSystemVector`. This includes copying the `DataCommunicator`, `DistributedNumbering`, local data, non-local data map, and the `DistributedVectorExporter` if it exists.

### Operators

`DistributedSystemVector` supports several operators. Most arithmetic operators work on the **local data** of each MPI process.

*   **Assignment (`=`):**
    ```cpp
    DistributedSystemVector& operator=(DistributedSystemVector const& rOtherVector);
    ```
    Assigns the local data from `rOtherVector` to the current vector. Requires `LocalSize()` to be the same for both vectors.

*   **Addition Assignment (`+=`):**
    ```cpp
    DistributedSystemVector& operator+=(const DistributedSystemVector& rOtherVector);
    ```
    Adds `rOtherVector`'s local data to the current vector's local data element-wise. Requires `LocalSize()` to be the same.

*   **Subtraction Assignment (`-=`):**
    ```cpp
    DistributedSystemVector& operator-=(const DistributedSystemVector& rOtherVector);
    ```
    Subtracts `rOtherVector`'s local data from the current vector's local data element-wise. Requires `LocalSize()` to be the same.

*   **Scalar Multiplication Assignment (`*=`):**
    ```cpp
    DistributedSystemVector& operator*=(const TDataType& multiplier_factor);
    ```
    Multiplies each element of the local data by `multiplier_factor`.

*   **Scalar Division Assignment (`/=`):**
    ```cpp
    DistributedSystemVector& operator/=(const TDataType& divide_factor);
    ```
    Divides each element of the local data by `divide_factor`.

*   **Element Access (`()` and `[]`):**
    ```cpp
    TDataType& operator()(IndexType I);
    const TDataType& operator()(IndexType I) const;
    TDataType& operator[](IndexType I);
    const TDataType& operator[](IndexType I) const;
    ```
    Provides read and write access to elements **within the local data partition** using their local index `I`. These operators do not access non-local data.

### Key Methods

Key methods for managing and operating on `DistributedSystemVector`:

*   **`GetComm() const`**
    ```cpp
    const DataCommunicator& GetComm() const;
    ```
    Returns a const reference to the `DataCommunicator` (e.g., `MPICommunicator`) used by this vector.

*   **`GetNumbering() const`**
    ```cpp
    const DistributedNumbering<IndexType>& GetNumbering() const;
    ```
    Returns a const reference to the `DistributedNumbering` object, which manages the mapping between global indices and local indices for each MPI process.

*   **`Clear()`**
    ```cpp
    void Clear();
    ```
    Clears the `mLocalData` vector, effectively setting its size to 0 on the current process.

*   **`SetValue(const TDataType value)`**
    ```cpp
    void SetValue(const TDataType value);
    ```
    Assigns the given `value` to all elements in the `mLocalData` (the local partition of the vector).

*   **`Size() const`**
    ```cpp
    IndexType Size() const;
    ```
    Returns the total global size of the distributed vector, as determined by the `DistributedNumbering`.

*   **`LocalSize() const`**
    ```cpp
    IndexType LocalSize() const;
    ```
    Returns the number of elements stored locally on the current MPI process.

*   **`GetLocalData()`**
    ```cpp
    DenseVector<TDataType>& GetLocalData();
    const DenseVector<TDataType>& GetLocalData() const;
    ```
    Returns a reference (or const reference) to the underlying `DenseVector<TDataType>` that stores the local portion of the vector's data.

*   **`GetExporter() const`**
    ```cpp
    const DistributedVectorExporter<TIndexType>& GetExporter() const;
    ```
    Returns a const reference to the `DistributedVectorExporter`. This object is responsible for managing the communication plan for non-local data during assembly (`FinalizeAssemble`). It is initialized in `BeginAssemble()`. An error is thrown if called before `BeginAssemble()`.

*   **`AddEntry(TIndexType GlobalI)`**
    ```cpp
    void AddEntry(TIndexType GlobalI); // WARNING: NOT THREADSAFE!
    ```
    Adds an entry for `GlobalI` to the `mNonLocalData` map if `GlobalI` is not local. This is used to inform the system that this process will contribute to this non-local degree of freedom.
    **Caution:** This method is not thread-safe and should be called before `BeginAssemble()`.

*   **`AddEntries(TIteratorType it_begin, TIteratorType it_end)`**
    ```cpp
    template<class TIteratorType>
    void AddEntries(TIteratorType it_begin, TIteratorType it_end); // WARNING: NOT THREADSAFE!
    ```
    Adds multiple entries to `mNonLocalData` based on the provided iterators.
    **Caution:** This method is not thread-safe and should be called before `BeginAssemble()`.

*   **`Norm() const`**
    ```cpp
    TDataType Norm() const;
    ```
    Calculates the L2 (Euclidean) norm of the entire distributed vector. This involves calculating the norm of the local data and then performing a sum-reduction across all MPI processes using `GetComm().SumAll()`.

*   **`Dot(const DistributedSystemVector& rOtherVector, MpiIndexType gather_on_rank=0)`**
    ```cpp
    TDataType Dot(const DistributedSystemVector& rOtherVector, MpiIndexType gather_on_rank=0);
    ```
    Computes the dot product of this vector with `rOtherVector`. The local dot product is computed, and then results are summed across MPI processes using `GetComm().Sum()`. The result is only guaranteed to be correct on `gather_on_rank`.

*   **`Add(const TDataType factor, const DistributedSystemVector& rOtherVector)`**
    ```cpp
    void Add(const TDataType factor, const DistributedSystemVector& rOtherVector);
    ```
    Performs the operation `this->mLocalData += factor * rOtherVector.mLocalData`. This is a purely local operation. Requires `LocalSize()` to be the same for both vectors.

*   **`BeginAssemble()`**
    ```cpp
    void BeginAssemble();
    ```
    Prepares the vector for assembly. Key actions include:
        *   Initializing the `DistributedVectorExporter` (`mpexporter`) based on the current contents of `mNonLocalData`. The structure of `mNonLocalData` is considered frozen after this call.
        *   Setting all values in `mNonLocalData` to zero, ready to accumulate incoming contributions.

*   **`FinalizeAssemble()`**
    ```cpp
    void FinalizeAssemble();
    ```
    Completes the assembly process. It uses the `mpexporter` (initialized in `BeginAssemble`) to exchange and sum the contributions accumulated in `mNonLocalData` with the owning processes. The local data (`mLocalData`) of the respective owning processes is updated with these summed non-local contributions.

*   **`Assemble(const TVectorType& rVectorInput, const TIndexVectorType& EquationId)`**
    ```cpp
    template<class TVectorType, class TIndexVectorType >
    void Assemble(
        const TVectorType& rVectorInput,
        const TIndexVectorType& EquationId
    );
    ```
    Assembles contributions from a local `rVectorInput` into the distributed vector.
        *   If an `EquationId[i]` is local (owned by the current process), `rVectorInput[i]` is atomically added to the corresponding entry in `mLocalData`.
        *   If an `EquationId[i]` is non-local, `rVectorInput[i]` is atomically added to the corresponding entry in `mNonLocalData`. These contributions are then communicated during `FinalizeAssemble()`.
    It's crucial that `mNonLocalData` has been populated (either via graph constructor or `AddEntry`/`AddEntries`) for all non-local `EquationId`s before `BeginAssemble` is called.

*   **`PrintInfo(std::ostream& rOStream) const` and `PrintData(std::ostream& rOStream) const`**
    ```cpp
    void PrintInfo(std::ostream& rOStream) const;
    void PrintData(std::ostream& rOStream) const;
    ```
    `PrintInfo` typically prints a header or general information. `PrintData` prints details about the vector on the current rank, including its local size, total size, the range of global IDs it owns, and the contents of `mLocalData`.

### Distributed Aspects

The `DistributedSystemVector` is fundamentally designed for parallel computing environments where data is partitioned across multiple MPI processes.

*   **Data Partitioning (Local vs. Non-Local):**
    *   Each MPI process "owns" a specific range of global indices. The data corresponding to these owned indices is stored in `mLocalData` (a `DenseVector`) on that process.
    *   When a process needs to contribute to a global index that it does not own (a non-local index), these contributions are temporarily stored in `mNonLocalData`. This is an `std::unordered_map<IndexType, TDataType>` where the key is the global index and the value is the accumulated contribution.
    *   The `DistributedNumbering` object is crucial for determining if a global index is local or non-local to the current process and for converting between global and local indices.

*   **Assembly Workflow:**
    The typical assembly process involves three main stages:
    1.  **`BeginAssemble()`:**
        *   This method prepares the vector for accumulating contributions.
        *   Crucially, it initializes the `DistributedVectorExporter` (`mpexporter`). The exporter analyzes the global indices present as keys in `mNonLocalData` and sets up a communication plan to send these contributions to their owning processes.
        *   The structure of `mNonLocalData` (i.e., which non-local global IDs are present) **must not change** after `BeginAssemble()` is called, as the communication plan depends on it.
        *   Values in `mNonLocalData` are typically set to zero at this stage to ensure fresh accumulation.
    2.  **`Assemble(const TVectorType& rVectorInput, const TIndexVectorType& EquationId)` (called multiple times):**
        *   For each entry in `rVectorInput` and its corresponding `EquationId`:
            *   If `EquationId[i]` is local to the current process (determined by `GetNumbering().IsLocal(EquationId[i])`), the value `rVectorInput[i]` is atomically added to `mLocalData[local_id]`.
            *   If `EquationId[i]` is non-local, the value `rVectorInput[i]` is atomically added to `mNonLocalData[EquationId[i]]`. A `KRATOS_DEBUG_ERROR_IF` check ensures that this `EquationId[i]` was already present in `mNonLocalData` (meaning it was anticipated, e.g., from graph construction or `AddEntry`).
    3.  **`FinalizeAssemble()`:**
        *   This method completes the assembly by performing the necessary MPI communication.
        *   The `DistributedVectorExporter` (`mpexporter`) uses the communication plan established in `BeginAssemble()` to send the accumulated values from each process's `mNonLocalData` to the processes that own those global indices.
        *   The receiving processes then add these contributions to their respective `mLocalData` entries.

*   **Role of `DistributedNumbering`:**
    *   Manages the partitioning of global indices across MPI processes.
    *   Provides methods like `IsLocal(GlobalId)`, `LocalId(GlobalId)`, and `GlobalId(LocalId)`.
    *   Defines the overall `Size()` of the vector and the `LocalSize()` for each process.
    *   Stores information about which process owns which range of global IDs (e.g., `GetCpuBounds()`).

*   **Role of `DistributedVectorExporter`:**
    *   Created in `BeginAssemble()` based on the non-local entries this process intends to contribute to.
    *   Manages the communication (MPI sends/receives) required in `FinalizeAssemble()` to transfer data from `mNonLocalData` of contributing processes to `mLocalData` of owning processes.
    *   The `DistributedVectorImporter` (used in the test file example) performs the reverse operation: gathering specified non-local data onto the current process.

### Code Examples

Below are examples demonstrating the usage of `DistributedSystemVector`. These examples assume an MPI environment is initialized.

#### **1. Construction**

```cpp
#include "containers/distributed_system_vector.h"
#include "containers/distributed_sparse_graph.h"
#include "containers/distributed_numbering.h"
#include "mpi/includes/mpi_data_communicator.h" // For DataCommunicator
#include "includes/parallel_environment.h"     // For ParallelEnvironment

// ... (Assume r_comm is an initialized DataCommunicator, e.g., from ParallelEnvironment)
const Kratos::DataCommunicator& r_comm = Kratos::ParallelEnvironment::GetDefaultDataCommunicator();
Kratos::DistributedSystemVector<>::IndexType local_size_per_process = 10;

// Construction using DistributedNumbering
// This creates a numbering where each process owns 'local_size_per_process' consecutive global IDs.
Kratos::DistributedNumbering<Kratos::DistributedSystemVector<>::IndexType> numbering(r_comm, local_size_per_process);
Kratos::DistributedSystemVector<> vec_from_numbering(numbering);
vec_from_numbering.SetValue(0.0); // Initialize local data

// Construction using DistributedSparseGraph (more common for FEM assembly)
// First, create and finalize a DistributedSparseGraph.
// 'dofs_bounds' would typically be determined by the mesh partitioning.
// Example: if total_dofs = 40, world_size = 2, my_rank = 0, dofs_bounds = {0, 20}
//          if total_dofs = 40, world_size = 2, my_rank = 1, dofs_bounds = {20, 40}
// This calculation usually happens inside utilities like ComputeBounds.
// For simplicity, let's assume 'dist_graph' is an already constructed and finalized DistributedSparseGraph.
Kratos::DistributedSparseGraph<Kratos::DistributedSystemVector<>::IndexType> dist_graph( ... );
Kratos::DistributedSystemVector<> vec_from_graph(dist_graph);
vec_from_graph.SetValue(1.0);
```

#### **2. Common Operations (Mostly Local)**

```cpp
#include "containers/distributed_system_vector.h"
#include "containers/distributed_numbering.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "includes/parallel_environment.h"
// ... (Assume r_comm and numbering are set up as above)
const Kratos::DataCommunicator& r_comm = Kratos::ParallelEnvironment::GetDefaultDataCommunicator();
Kratos::DistributedSystemVector<>::IndexType local_size = 4;
Kratos::DistributedNumbering<Kratos::DistributedSystemVector<>::IndexType> numbering(r_comm, local_size);
Kratos::DistributedSystemVector<> a(numbering);
a.SetValue(5.0); // Each process sets its local elements to 5.0
Kratos::DistributedSystemVector<> b(numbering);
b.SetValue(3.0); // Each process sets its local elements to 3.0
// Element access (local indexing)
if (a.LocalSize() > 0) {
    a[0] = 10.0; // Modifies the first local element on each process
}
// Arithmetic operations (act on local data)
Kratos::DistributedSystemVector<> c(a); // Copy constructor
c += b; // c's local data becomes a's local data + b's local data (e.g., 10+3, 5+3, ...)
c -= b; // c's local data back to a's original state (e.g., 10, 5, ...)
c.Add(2.0, b); // c_local = c_local + 2.0 * b_local (e.g., 10 + 2*3, 5 + 2*3, ...)
c *= 0.5;
c /= 2.0;
// Global operations (involve communication)
double norm_c = c.Norm(); // Calculates global L2 norm
double dot_ab = a.Dot(b); // Calculates global dot product (result on rank 0 by default)
// Access local data directly
Kratos::DenseVector<double>& local_vector_a = a.GetLocalData();
if (local_vector_a.size() > 1) {
    local_vector_a[1] = 7.0; // Directly modifies the second local element of 'a'
}
// std::cout << "Vector a on rank " << r_comm.Rank() << ": ";
// a.PrintData(std::cout); // Prints local data for each rank
```

#### **3. Distributed Assembly Process**

This example is based on `DistributedSystemVectorAssembly` test.

```cpp
#include "containers/distributed_system_vector.h"
#include "containers/distributed_sparse_graph.h"
#include "mpi/testing/mpi_testing.h" // For test utilities if needed
#include "mpi/tests/test_utilities/distributed_sparse_containers_test_utilities.h" // For ComputeBounds, ElementConnectivities

// ...

const auto& r_comm = Kratos::ParallelEnvironment::GetDefaultDataCommunicator();
int world_size = r_comm.Size();
int my_rank = r_comm.Rank();

// 1. Define DOF distribution (partitioning)
// Example: 40 total DOFs, distributed among 'world_size' processes
auto dofs_bounds = Kratos::Testing::DistributedSparseContainersTestUtilities::ComputeBounds<Kratos::DistributedSystemVector<>::IndexType>(40, world_size, my_rank);

// 2. Define element connectivities (local to each process)
// Example: 31 total elements, 'el_bounds' defines which elements this process handles
auto el_bounds = Kratos::Testing::DistributedSparseContainersTestUtilities::ComputeBounds<Kratos::DistributedSystemVector<>::IndexType>(31, world_size, my_rank);
const auto connectivities = Kratos::Testing::DistributedSparseContainersTestUtilities::ElementConnectivities(el_bounds);
// 'connectivities' is typically std::vector<std::vector<GlobalDofIdType>>

// 3. Create and Finalize DistributedSparseGraph
// The graph needs to know about all local and non-local DOFs this process will interact with.
Kratos::DistributedSparseGraph<Kratos::DistributedSystemVector<>::IndexType> dist_sparse_graph(dofs_bounds[1] - dofs_bounds[0], r_comm);
Kratos::IndexPartition<Kratos::DistributedSystemVector<>::IndexType>(connectivities.size()).for_each([&](Kratos::DistributedSystemVector<>::IndexType i) {
    dist_sparse_graph.AddEntries(connectivities[i]);
});
dist_sparse_graph.Finalize(); // Finalizes graph structure, prepares communication patterns

// 4. Create DistributedSystemVector using the graph
Kratos::DistributedSystemVector<> b_rhs(dist_sparse_graph);
b_rhs.SetValue(0.0); // Initialize vector (local parts) to zero

// 5. Assembly Loop
b_rhs.BeginAssemble(); // Prepares for assembly, initializes exporter based on mNonLocalData (from graph)

for(const auto& elem_conn : connectivities) { // Loop over elements local to this process
    // In a real FEM code, 'elem_contribution_vector' would be computed by an Element.
    // For this example, assume each DOF in the element gets a +1.0 contribution.
    Kratos::Vector elem_contribution_vector(elem_conn.size());
    for(size_t i=0; i<elem_conn.size(); ++i) {
        elem_contribution_vector[i] = 1.0;
    }

    // Assemble local contributions into the global vector 'b_rhs'
    // 'elem_conn' contains the global DOF IDs for this element.
    b_rhs.Assemble(elem_contribution_vector, elem_conn);
}

b_rhs.FinalizeAssemble(); // Exchanges and sums non-local contributions via MPI

// Now b_rhs contains the assembled vector, distributed across processes.
// Verification (example from test):
// Compare local values with a reference map.
// auto reference_b_map = Kratos::Testing::DistributedSparseContainersTestUtilities::GetReferenceVectorAsMap(dofs_bounds);
for(unsigned int i = 0; i < b_rhs.LocalSize(); ++i) {
    auto global_i = b_rhs.GetNumbering().GlobalId(i);
    auto it = reference_b_map.find(global_i);
    KRATOS_EXPECT_NEAR(b_rhs(i),  it->second, 1e-14 );
}
```