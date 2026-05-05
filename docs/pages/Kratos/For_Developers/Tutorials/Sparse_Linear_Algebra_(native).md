---
title: Sparse Linear Algebra (native Kratos implementation)
keywords:
tags: [Sparse_Linear_Algebra_(native).md]
sidebar: kratos_for_developers
summary:
---

# **Sparse Linear Algebra (native Kratos implementation)**

The Kratos platform provides a basic native implementation of sparse linear algebra tools. These tools are primarily designed to facilitate Finite Element Method (FEM) operations and to offer a consistent interface for both serial (Shared Memory Parallelism, SMP) and distributed memory (Message Passing Interface, MPI) parallel computing environments.

## Key Data Structures

For detailed information on the specific data structures discussed in this tutorial, please refer to their individual documentation pages:

### SMP
*   [`CsrMatrix`](../Data_Structures/CsrMatrix.md)
*   [`SparseGraph`](../Data_Structures/SparseGraph.md)
*   [`SystemVector`](../Data_Structures/SystemVector.md)

### MPI
*   [`DistributedCsrMatrix`](../Data_Structures/DistributedCsrMatrix.md)
*   [`DistributedSparseGraph`](../Data_Structures/DistributedSparseGraph.md)
*   [`DistributedSystemVector`](../Data_Structures/DistributedSystemVector.md)

## Introduction to Kratos Sparse Linear Algebra
The Kratos sparse linear algebra library offers fundamental capabilities for constructing and manipulating sparse matrices, particularly focusing on the Compressed Sparse Row (CSR) format, which is commonly used in FEM applications. A key aspect of CSR matrix construction is defining the matrix's sparsity pattern (i.e., which entries are non-zero). In Kratos, this is typically achieved using "graph" objects. These graph objects allow for the dynamic or static definition of non-zero entries before the matrix itself is numerically populated.

## Serial Sparse Linear Algebra (Shared Memory Parallelism, SMP)

### Graph Objects in SMP ([`SparseGraph`](../Data_Structures/SparseGraph.md), `SparseContiguousRowGraph`)
In the SMP (Shared Memory Parallelism, typically leveraging multi-core processors via threads like OpenMP) environment, Kratos offers two primary "graph-type" objects to define the sparsity pattern of matrices. These objects are crucial for efficient CSR matrix construction.

*   **[`SparseGraph`](../Data_Structures/SparseGraph.md)**: This object is typically used when the matrix structure is dynamic or the total number of rows is not known beforehand. It allows for flexible addition of entries.
    *   It allows constructing a graph where neither the number of rows nor the number of columns needs to be pre-determined.
    *   Entries are added using methods like `AddEntry(I,J)` for a single entry, or `AddEntries(list_of_ids)` to add all combinations for a given list of IDs (creating a dense block), or `AddEntries(row_ids, column_ids)` for rectangular blocks.
    *   After all entries are added, `Finalize()` must be called. This prepares the graph for use in constructing a `CsrMatrix`, typically by sorting and removing duplicate entries.
    *   Crucially, [`SparseGraph`](../Data_Structures/SparseGraph.md) is **not threadsafe**. If different parts of the graph are constructed in parallel, they must be merged subsequently.

*   **`SparseContiguousRowGraph`**: This object is preferred when the matrix structure is fixed, meaning the number of rows is known and will not change.
    *   It requires the number of rows to be specified at construction.
    *   Similar to `SparseGraph`, it uses `AddEntry` and `AddEntries` methods for populating the graph. The `Finalize()` call is also necessary.
    *   A key advantage is that `SparseContiguousRowGraph` **is threadsafe** for the assembly process (i.e., calls to `AddEntries`). This makes it suitable for parallel construction of the graph where different threads might add entries to different rows, or even the same row, concurrently.

The general workflow for graph construction involves creating an instance of either graph type and then iteratively adding entries based on element connectivities (i.e., the set of global Degree of Freedom (DOF) indices associated with each element) or other criteria.

The following C++ snippet illustrates a typical loop for populating a graph (using `SparseGraph` as an example, but `SparseContiguousRowGraph` would be similar if the number of rows were known at `Agraph`'s construction):
~~~cpp
    [`SparseGraph`](../Data_Structures/SparseGraph.md)<> Agraph; // Alternatively, use `SparseContiguousRowGraph` if num_rows is known
    for(const auto& c : connectivities) // c is a list of global DOF IDs, e.g., from an element
        Agraph.AddEntries(c); // Adds all (i,j) for i,j in c. For c=[1,2], adds (1,1),(1,2),(2,1),(2,2)
    Agraph.Finalize(); // Marks the graph as complete and ready for CsrMatrix construction
~~~

Both graph types offer these common functionalities (where `I` and `J` refer to global indices):
*   **AddEntry(I,J)**: Adds a single entry (I,J) to the graph.
*   **AddEntries(list_of_ids)**: Adds all entries (list_of_ids(i), list_of_ids(j)) for all global indices i and j within the provided list.
*   **AddEntries(row_ids, column_ids)**: Adds entries for a rectangular block defined by lists of global row and column indices.
*   **Has(I,J)**: Checks if a specific entry (I,J) is present in the graph after finalization.
*   **Size()**: Returns the number of rows in the graph.

Graphs can be iterated to access the non-zero positions (I,J) as shown below:
~~~cpp
    for(auto it=rAgraph.begin(); it!=rAgraph.end(); ++it)
    {
        const auto I = it.GetRowIndex(); // Gets the row index
        for(auto J : *it ) // Iterates over column indices in that row
        {
            // Access entry (I,J)
        }
    }
~~~

Additionally, these graph objects provide utility functions:
*   **ExportCSRArrays**: Converts the graph representation into raw CSR arrays (row pointers and column indices). This is useful for interfacing with other libraries or for direct `CsrMatrix` construction.
*   **ExportSingleVectorRepresentation**: Provides a compact, single-vector representation of the graph, useful for serialization or fast graph reconstruction.

### [`CsrMatrix`](../Data_Structures/CsrMatrix.md)
The Kratos native implementation of a CSR (Compressed Sparse Row) matrix, [`CsrMatrix`](../Data_Structures/CsrMatrix.md), is optimized for typical Finite Element Method (FEM) assembly operations. Once a graph defining the sparsity pattern (e.g., from a `SparseGraph` or `SparseContiguousRowGraph`) is finalized, a `CsrMatrix` can be constructed from it.

**Construction from a Graph:**
The primary way to create a `CsrMatrix` is by passing a finalized graph object to its constructor. This graph dictates the non-zero structure of the matrix.
```cpp
    // Assuming Agraph is a finalized SparseGraph or SparseContiguousRowGraph
    CsrMatrix<double> A(Agraph);
```
This constructor allocates the necessary memory for the matrix based on the graph's structure. The matrix values are typically initialized to zero.

**Assembly Workflow:**
Populating the `CsrMatrix` with values follows a specific workflow, designed to be efficient and threadsafe for parallel assembly.

1.  **`BeginAssemble()`**: This method must be called before any assembly operations begin. It prepares the matrix for assembly, for instance, by initializing values to zero or setting up internal data structures.
2.  **`Assemble()`**: This is the core method for adding values to the matrix. It has several overloads:
    *   `Assemble(value, I, J)`: Assembles a single scalar `value` into the matrix at position (I,J). This position must exist in the graph from which the matrix was constructed.
    *   `Assemble(local_matrix, positions)`: This is the standard FEM assembly. It assembles a (typically dense) `local_matrix` (e.g., an element stiffness matrix) into the global matrix according to the `positions` (a list of global indices, e.g., equation IDs or DOFs, associated with the element).
    *   `Assemble(rectangular_matrix, row_ids, col_ids)`: Assembles a `rectangular_matrix` into the global matrix using specified lists of global `row_ids` and `col_ids`.
    The assembly operations (calls to `Assemble`) on **[`CsrMatrix`](../Data_Structures/CsrMatrix.md)** objects are **threadsafe**. This means multiple threads can call `Assemble()` concurrently on the same `CsrMatrix` object, which is crucial for parallelizing the assembly phase of FEM computations.
3.  **`FinalizeAssemble()`**: After all contributions have been assembled, this method must be called. It finalizes the assembly process, which might involve summing up contributions from different threads if they were writing to intermediate storage, or other necessary cleanup operations.

A typical assembly loop looks like this:
```cpp
    CsrMatrix<double> A(Agraph); // Construct from an existing graph

    A.BeginAssemble(); // Prepare for assembly
    // This loop can often be parallelized (e.g., using OpenMP)
    for(const auto& element_data : elements_list)
    {
        LocalMatrixType local_stiffness_matrix; // Calculated for the current element
        EquationIdVectorType global_dof_ids; // Global DOF IDs for the current element
        // ... calculate local_stiffness_matrix and global_dof_ids ...
        A.Assemble(local_stiffness_matrix, global_dof_ids);
    }
    A.FinalizeAssemble(); // Finalize assembly
```

The class also implements these 3 entry points for assembly, as mentioned (where `I,J` and IDs in `positions` are global indices):
* **Assemble(value,I,J)**: Assembling the single scalar "value" into I,J (I,J needs to be in the graph).
* **Assemble(square_matrix,positions)**: "Standard" FEM assembly of a square_matrix using global indices specified in `positions`.
* **Assemble(rectangular_matrix,row_ids,col_ids)**: Assembling of a rectangular matrix using global `row_ids` and `col_ids`.


**Key Methods and Usage:**
Beyond assembly, `CsrMatrix` provides several useful methods:

*   **`SpMV(x, y)`**: Performs the sparse matrix-vector product `y += A*x`. This is a fundamental operation in iterative solvers.
*   **`SpMV(alpha, x, beta, y)`**: Performs a generalized sparse matrix-vector product `y = alpha*A*x + beta*y`.
*   **`TransposeSpMV(x, y)`**: Performs the transpose sparse matrix-vector product `y += A^T*x`.
*   **`TransposeSpMV(alpha, x, beta, y)`**: Performs the generalized transpose sparse matrix-vector product `y = alpha*A^T*x + beta*y`.
*   **Data Accessors**:
    *   `index1_data()`: Returns a pointer/reference to the underlying array of row pointers (CSR `row_ptr` array).
    *   `index2_data()`: Returns a pointer/reference to the underlying array of column indices (CSR `col_ind` array).
    *   `value_data()`: Returns a pointer/reference to the underlying array of non-zero values.
    These accessors are useful for advanced operations, interfacing with external libraries, or for direct manipulation of the matrix data, but should be used with caution as they bypass the class's encapsulation.
*   **`size1()`**: Returns the number of rows.
*   **`size2()`**: Returns the number of columns.
*   **`nnz()`**: Returns the number of non-zero entries.
*   **`ToMap()`**: Can be useful for debugging or inspection, returning a map-like representation `{(I,J): value}` of the sparse matrix entries. (Note: This can be memory intensive for large matrices).

### System Vectors in SMP ([`SystemVector`](../Data_Structures/SystemVector.md))
In the context of serial (SMP) computations, the [`SystemVector`](../Data_Structures/SystemVector.md) is the standard Kratos data structure used to represent vectors, such as the Right-Hand Side (RHS) or solution vectors in a system of linear equations `Ax=b`. It is essentially a wrapper around a standard C++ vector, extended with functionalities for FEM assembly and common vector operations.

**Construction and Sizing:**
A `SystemVector` can be constructed in several ways:
*   **Default constructor**: Creates an empty vector.
    ```cpp
    SystemVector<double> b;
    ```
*   **Sized constructor**: Creates a vector of a specified size, typically initialized to zeros. This is common when the size is known, for example, from the number of rows in a matrix or a graph.
    ```cpp
    unsigned int vector_size = A.size1(); // A is a CsrMatrix or a graph object
    SystemVector<double> x(vector_size);
    ```
*   **Copy constructor/Assignment**: Can be created from another `SystemVector` or compatible Kratos/Trilinos vector types.

**Assembly Workflow:**
Similar to `CsrMatrix`, `SystemVector` has an assembly workflow, which is crucial for constructing vectors like the RHS in FEM applications. This workflow is also designed to be threadsafe.

1.  **`BeginAssemble()`**: This call prepares the vector for assembly. It typically initializes all entries to zero. This is essential to ensure that contributions from different elements or sources are summed correctly without interference from previous values.
2.  **`Assemble()` / `AssembleEntry()`**:
    *   `Assemble(local_vector, positions)`: This is the most common method for FEM. It adds the entries of a `local_vector` (e.g., an element load vector) to the global `SystemVector` at the global indices specified by `positions` (e.g., global equation IDs corresponding to the element's DOFs).
    *   `AssembleEntry(value, I)`: Adds a single `value` to the entry at global index `I` of the `SystemVector`. This is useful for applying point loads or other single-entry modifications.
    The assembly operations on `SystemVector` are **threadsafe**, allowing multiple threads to contribute to the vector concurrently.
3.  **`FinalizeAssemble()`**: This call finalizes the assembly process. While often a no-op in the serial `SystemVector` (as direct summation is usually performed), it's good practice to include it for consistency with the distributed `DistributedSystemVector` and to ensure any potential future finalization steps are handled.

Example of assembling an RHS vector:
```cpp
    SystemVector<double> b(Agraph.Size()); // Size based on graph's number of rows

    b.BeginAssemble();
    // This loop can be parallelized
    for(const auto& element_data : elements_list) // Assuming elements_list holds element data
    {
        LocalVectorType element_rhs_contribution; // Calculated for the current element
        EquationIdVectorType global_dof_ids;      // Global DOF IDs for the current element
        // ... calculate element_rhs_contribution and global_dof_ids ...
        b.Assemble(element_rhs_contribution, global_dof_ids);
    }
    // Example of adding a single entry
    IndexType specific_global_dof_id = ...; // A specific global DOF ID
    double load_value = 10.0;
    b.AssembleEntry(load_value, specific_global_dof_id);

    b.FinalizeAssemble();
```

**Other Useful Operations:**
`SystemVector` supports a range of common vector operations:
*   **`Norm()` / `Norm2()`**: Calculates the Euclidean (L2) norm of the vector.
*   **`Dot(other_vector)`**: Computes the dot product with another `SystemVector`.
*   **`SetValue(value)`**: Sets all entries of the vector to a specified `value`.
*   **Element Access**: Standard C++ `operator[]` can be used to access or modify individual entries, e.g., `x[i] = value;` or `value = x[i];`.
*   Standard arithmetic operations like addition, subtraction, and scaling are also available, often through operator overloads or dedicated methods.

The `SystemVector` provides the necessary functionality for handling global vectors in serial FEM computations, from assembly to their use in linear solvers and post-processing.

### SMP Assembly Workflow Example
This section provides a concise C++ example demonstrating a typical serial (SMP) workflow using the previously discussed `SparseContiguousRowGraph`, `CsrMatrix`, and `SystemVector` classes. This example simulates a very simple Finite Element problem with a few elements and degrees of freedom (DOFs).

```cpp
// Assuming necessary Kratos headers are included:
#include "includes/kratos_parameters.h"
#include "sparse_contiguous_row_graph.h" // Or "sparse_graph.h"
#include "csr_matrix.h"
#include "system_vector.h"
#include <vector>

// Define some Kratos types for clarity (actual types might vary)
using IndexType = std::size_t;
using ValueType = double;
using LocalMatrixType = Kratos::Matrix; // Assuming a Kratos dense matrix for local systems
using LocalVectorType = Kratos::Vector; // Assuming a Kratos dense vector for local systems
using EquationIdVectorType = std::vector<IndexType>;

void SmpAssemblyExample()
{
    // 1. Define Problem Size (e.g., number of DOFs)
    const IndexType number_of_dofs = 5; // Small example: 5 DOFs

    // 2. Create a Graph (using SparseContiguousRowGraph for fixed size and thread-safety)
    // For dynamic graphs or when size is unknown initially, SparseGraph could be used.
    Kratos::SparseContiguousRowGraph<IndexType> fem_graph(number_of_dofs);

    // 3. Populate the Graph with sample connectivities
    // These represent the DOFs involved in each "element".
    // Element 1: DOFs 0, 1, 2
    EquationIdVectorType elem1_dofs = {0, 1, 2};
    fem_graph.AddEntries(elem1_dofs);

    // Element 2: DOFs 1, 2, 3
    EquationIdVectorType elem2_dofs = {1, 2, 3};
    fem_graph.AddEntries(elem2_dofs);

    // Element 3: DOFs 2, 3, 4
    EquationIdVectorType elem3_dofs = {2, 3, 4};
    fem_graph.AddEntries(elem3_dofs);

    // 4. Finalize the Graph
    // This prepares the graph for CsrMatrix construction (sorts, removes duplicates, etc.)
    fem_graph.Finalize();

    // 5. Create a CsrMatrix from the finalized graph
    Kratos::CsrMatrix<ValueType, IndexType> A(fem_graph);

    // 6. Create SystemVectors for RHS (b) and solution (x)
    Kratos::SystemVector<ValueType, IndexType> b(A.size1()); // RHS vector, size based on matrix rows
    Kratos::SystemVector<ValueType, IndexType> x(A.size1()); // Solution vector, initialized to zeros usually

    // Initialize x with some dummy values for the SpMV example
    for(IndexType i = 0; i < x.size(); ++i) {
        x[i] = static_cast<ValueType>(i + 1);
    }

    // 7. Assemble Global Matrix (A) and RHS Vector (b)
    A.BeginAssemble();
    b.BeginAssemble();

    // --- Loop over conceptual "elements" ---
    // Element 1 (DOFs: 0, 1, 2)
    LocalMatrixType K_elem1(3, 3); // 3x3 local stiffness matrix
    LocalVectorType f_elem1(3);    // 3x1 local force vector
    // Fill K_elem1 and f_elem1 with some dummy values
    for(int i=0; i<3; ++i) { for(int j=0; j<3; ++j) K_elem1(i,j) = (i==j) ? 2.0 : -1.0; } // Simple Laplacian-like
    for(int i=0; i<3; ++i) { f_elem1(i) = 1.0 * (i+1); }

    A.Assemble(K_elem1, elem1_dofs);
    b.Assemble(f_elem1, elem1_dofs);

    // Element 2 (DOFs: 1, 2, 3)
    LocalMatrixType K_elem2(3, 3);
    LocalVectorType f_elem2(3);
    for(int i=0; i<3; ++i) { for(int j=0; j<3; ++j) K_elem2(i,j) = (i==j) ? 2.0 : -0.5; }
    for(int i=0; i<3; ++i) { f_elem2(i) = 0.5 * (i+1); }

    A.Assemble(K_elem2, elem2_dofs);
    b.Assemble(f_elem2, elem2_dofs);

    // Element 3 (DOFs: 2, 3, 4)
    LocalMatrixType K_elem3(3, 3);
    LocalVectorType f_elem3(3);
    for(int i=0; i<3; ++i) { for(int j=0; j<3; ++j) K_elem3(i,j) = (i==j) ? 3.0 : -1.2; }
    for(int i=0; i<3; ++i) { f_elem3(i) = 1.2 * (i+1); }

    A.Assemble(K_elem3, elem3_dofs);
    b.Assemble(f_elem3, elem3_dofs);
    // --- End of element loop ---

    A.FinalizeAssemble();
    b.FinalizeAssemble();

    // 8. Perform a Sparse Matrix-Vector Product (e.g., y = A*x)
    // This is a common operation in solvers. Here, we use b as the result y (b = A*x).
    // For a real solve, b would be the RHS and x the unknown solution.
    // Let's re-initialize b for this operation to make it clear:
    b.SetValue(0.0); // Set all entries of b to 0 before y = A*x type of operation
                     // or use SpMV(const VectorType& rX, VectorType& rY, bool TransposeA = false)
                     // if you want y = A*x
                     // The existing SpMV(x,b) does b += A*x.
                     // For b = A*x, if b is not zero, first set b=0 or use a different SpMV overload if available.
    // For this example, let's compute y_result = A * x.
    // If SpMV is y += A*x, y_result must be zeroed first.
    Kratos::SystemVector<ValueType, IndexType> y_result(A.size1());
    y_result.SetValue(0.0); // Ensure y_result is zero before y_result += A*x
    A.SpMV(x, y_result);

    // For demonstration, let's print out the resulting vector y_result
    // Note: In a real application, use Kratos' logger (KRATOS_WATCH) for output.
    std::cout << "Resulting vector y_result after SpMV: " << y_result << std::endl;
}
```
This example demonstrates the end-to-end process:
1.  **Graph Construction**: A `SparseContiguousRowGraph` is created and populated with Degree of Freedom (DOF) connectivities from three simple elements. This defines the non-zero pattern of the matrix.
2.  **Matrix and Vector Creation**: A `CsrMatrix` `A` is created from the graph. `SystemVector`s `b` (for RHS) and `x` (for solution) are also created.
3.  **Assembly**: The code simulates looping through elements, calculating their local stiffness matrices (`K_elem`) and local force vectors (`f_elem`), and then assembling them into the global matrix `A` and vector `b` using their respective `Assemble` methods with global DOF IDs.
4.  **Matrix Operation**: Finally, it shows a conceptual `SpMV` (Sparse Matrix-Vector product) operation, `y_result = A*x`. In a real solver, `x` would typically be the current solution guess and `b` the right-hand side, or `A*x` would be used in residual calculations.

This workflow is fundamental to many FEM-based simulations where system matrices and vectors are built by accumulating contributions from individual elements. The use of `SparseContiguousRowGraph` makes graph population thread-safe, and the assembly methods of `CsrMatrix` and `SystemVector` are also designed for thread-safety in SMP environments.

## Distributed Sparse Linear Algebra (Message Passing Interface, MPI)

### Core MPI Concepts (`DataCommunicator`, `DistributedNumbering`)
Distributed sparse linear algebra in Kratos heavily relies on the Message Passing Interface (MPI) for inter-process communication. Two fundamental concepts underpin this: the `DataCommunicator` and `DistributedNumbering` (often represented by `RowNumbering` and `ColNumbering` objects internally within distributed classes).

#### `DataCommunicator`
The `DataCommunicator` is Kratos's abstraction layer over MPI. It handles all MPI communication calls (e.g., send, receive, broadcast, gather, scatter, reductions) between different processes (MPI ranks).
*   **Role**: It provides a consistent interface for parallel communication, regardless of whether Kratos is compiled with MPI support or running in serial mode (where it acts as a no-op for most MPI calls). This allows writing parallel-agnostic code to a large extent.
*   **Obtaining a `DataCommunicator`**: Typically, you get access to a `DataCommunicator` object that represents the global MPI communicator (e.g., `MPI_COMM_WORLD`). In Kratos, this is often retrieved via the `ParallelEnvironment`:
    ```cpp
    // Get the default DataCommunicator for the current parallel environment
    const Kratos::DataCommunicator& rComm = Kratos::ParallelEnvironment::GetDefaultDataCommunicator();
    // Or from a ModelPart in an MPI simulation:
    // const Kratos::DataCommunicator& rComm = rModelPart.GetCommunicator();
    ```
    This `rComm` object is then passed to distributed data structures like `DistributedSparseGraph`, `DistributedCsrMatrix`, and `DistributedSystemVector` to enable them to perform necessary communications.

#### `DistributedNumbering` (`RowNumbering`, `ColNumbering`)
The concept of `DistributedNumbering` is crucial for managing data distribution and mapping between global and local views in a parallel environment. While there isn't a single class named `DistributedNumbering`, this role is typically fulfilled by objects (often internally within distributed matrices/vectors, or as standalone utilities if needed) that manage how global indices (e.g., Degrees of Freedom (DOFs), matrix rows/columns) are partitioned and accessed across MPI ranks. For `DistributedCsrMatrix`, these are often referred to as `RowNumbering` and `ColNumbering`.

*   **Purpose**:
    *   **Data Ownership**: Defines which MPI rank "owns" which global Degree of Freedom (DOF) or matrix/vector row. Each global ID is uniquely owned by one rank.
    *   **Global-to-Local Mapping**: Translates a global ID (unique across all ranks) to a local ID (unique within a specific rank, typically for accessing data in local arrays).
    *   **Ghost/Non-Local Data**: Identifies "ghost" or "non-local" IDs. These are IDs that a rank needs to know about (e.g., for assembling contributions from shared interface nodes) but does not own. The numbering system helps manage communication for these ghost entities.

*   **Construction**:
    A distributed numbering is typically established by providing information about how many global IDs each rank owns. For example, one common way is to define the size of the partition (number of rows/DOFs) handled by each rank.
    ```cpp
    // Example: Building a row numbering (conceptual)
    // This information is often derived from the partitioning of a ModelPart
    IndexType number_of_local_rows_on_my_rank = ...;
    const DataCommunicator& rComm = ParallelEnvironment::GetDefaultDataCommunicator();

    // Internally, distributed objects use this info to build their numbering.
    // For instance, when constructing a DistributedSparseGraph:
    // DistributedSparseGraph<IndexType> dist_graph(number_of_local_rows_on_my_rank, rComm);
    // This dist_graph will then internally create its RowNumbering.
    ```
    The `ColNumbering` for a `DistributedCsrMatrix` is often derived from the `RowNumbering` and the graph structure, identifying which column indices are local (part of the diagonal block) and which are remote (part of the off-diagonal block).

*   **Key Functionalities**:
    Objects managing distributed numbering provide key functionalities:
    *   `IsLocal(GlobalId)`: Returns `true` if the given `GlobalId` is owned by the current MPI rank.
    *   `OwnerRank(GlobalId)`: Returns the MPI rank that owns the `GlobalId`. For a local ID, this would be the current rank.
    *   `GlobalToLocal(GlobalId)`: Converts a `GlobalId` (that is local to the current rank) to its corresponding `LocalId`.
    *   `LocalToGlobal(LocalId)`: Converts a `LocalId` on the current rank back to its `GlobalId`.
    *   Methods to get lists of local IDs, ghost IDs, etc.

Understanding `DataCommunicator` and the principles of `DistributedNumbering` (how data is split and referred to across ranks) is essential for working with distributed matrices and vectors in Kratos, as these components handle the complexities of parallel data management and communication. The `DistributedSparseGraph`, `DistributedCsrMatrix`, and `DistributedSystemVector` classes encapsulate much of this logic.

### Distributed Graph Objects ([`DistributedSparseGraph`](../Data_Structures/DistributedSparseGraph.md))
The [`DistributedSparseGraph`](../Data_Structures/DistributedSparseGraph.md) is the cornerstone for defining the sparsity pattern of distributed matrices in an MPI environment. It handles the complexities of distributed data management across multiple MPI ranks.

**Construction:**
A `DistributedSparseGraph` is constructed by providing:
1.  The number of **locally owned rows** that this MPI rank will manage. This implicitly defines a part of the `RowNumbering`.
2.  A `DataCommunicator` object, which the graph will use for all necessary MPI communications.

The constructor typically takes the local size (number of rows owned by the current rank) and the `DataCommunicator`. The global numbering of rows is established collectively based on the local sizes provided by each rank.
```cpp
    // Get the DataCommunicator
    const DataCommunicator& rComm = ParallelEnvironment::GetDefaultDataCommunicator();

    // Determine the number of rows locally owned by this rank
    IndexType local_number_of_rows = ...; // This usually comes from the model part partition

    // Construct the DistributedSparseGraph
    DistributedSparseGraph<IndexType> dist_graph(local_number_of_rows, rComm);
```
Internally, this constructor sets up the `RowNumbering` which determines which global row indices are owned by which rank and how they map to local indices.

**Adding Entries (`AddEntry`, `AddEntries`):**
Adding entries to a `DistributedSparseGraph` is similar to its serial counterpart, but with added logic for distribution:
*   All entries are specified using **global indices**.
*   **Locally Owned Rows**: If an entry `(I,J)` is added where row `I` is owned by the current MPI rank, the entry `(I,J)` is stored directly in the local graph structure of this rank.
*   **Non-Locally Owned Rows (Remote Rows)**: If an entry `(I,J)` is added where row `I` is owned by a *different* MPI rank, this entry is temporarily cached by the current rank. These cached entries are destined for other ranks.
*   The `AddEntry(I,J)`, `AddEntries(list_of_ids)`, and `AddEntries(row_ids, column_ids)` methods all operate with global indices and handle this local storage vs. remote caching internally.

This process allows each rank to build its portion of the graph based on its local elements or conditions, without needing to know upfront which other ranks might also contribute to its owned rows or require entries from its rows.
```cpp
    // Global connectivities for an element processed by this rank
    EquationIdVectorType global_connectivities = {global_dof_id_1, global_dof_id_2, global_dof_id_3};

    // Add entries to the graph using global IDs
    // The graph internally determines if these correspond to local or remote rows.
    dist_graph.AddEntries(global_connectivities);
```

**Finalization (`Finalize()`):**
The `Finalize()` method is a **collective operation** and is of critical importance for `DistributedSparseGraph`.
*   When `Finalize()` is called, all MPI ranks communicate to exchange the cached non-local entries. Each rank sends the entries it has cached for other ranks, and receives entries that other ranks have cached for its owned rows.
*   After this communication, each rank integrates the received entries into its local graph structure.
*   The method also typically sorts and removes duplicate entries within each rank's local graph portion.
*   **All MPI ranks must call `Finalize()` concurrently.** It acts as a synchronization point, ensuring that the complete distributed graph structure is consistent across all processes before proceeding (e.g., to construct a `DistributedCsrMatrix`).
```cpp
    // ... after all AddEntry/AddEntries calls on all ranks ...
    dist_graph.Finalize(); // Collective MPI communication and graph finalization happens here
```

**Accessing the Local Graph:**
Once finalized, each rank holds the portion of the graph corresponding to its locally owned rows. This local part can be accessed for inspection or iteration if needed:
```cpp
    const auto& local_graph = dist_graph.GetLocalGraph();
    for(auto it = local_graph.begin(); it != local_graph.end(); ++it)
    {
        const auto I_local = it.GetRowIndex(); // This is a local row index
        // To get the global row index: dist_graph.RowNumbering().LocalToGlobal(I_local);
        for(auto J_global : *it ) // Column indices are typically stored as global
        {
            // Process entry (GlobalID(I_local), J_global)
        }
    }
```
The method `Size()` on a `DistributedSparseGraph` returns the total number of rows across all ranks (the global size), while iteration as shown above is typically over the local rows. The method `Has(I,J)` also allows querying for global entries. The `AddEntry`, `AddEntries` methods are available as in the serial case.

The `DistributedSparseGraph` thus provides the necessary distributed data structure to define the sparsity pattern for large-scale parallel computations, abstracting many of the MPI communication details from the user.

### Distributed CSR Matrices ([`DistributedCsrMatrix`](../Data_Structures/DistributedCsrMatrix.md))
The [`DistributedCsrMatrix`](../Data_Structures/DistributedCsrMatrix.md) class represents a CSR matrix partitioned across multiple MPI ranks. Each rank stores only its locally owned rows, but these rows can contain entries whose column indices point to DOFs owned by other ranks.

**Construction:**
A `DistributedCsrMatrix` is primarily constructed from a finalized [`DistributedSparseGraph`](../Data_Structures/DistributedSparseGraph.md). The graph provides the complete sparsity pattern (both local and remote connections) and the necessary `DataCommunicator` and `DistributedNumbering` (row numbering) information.
```cpp
    // Assuming 'dist_graph' is a finalized DistributedSparseGraph
    DistributedCsrMatrix<ValueType, IndexType> A(dist_graph);
    // The DataCommunicator and row numbering are inherited from the graph.
    // MPI communications may occur here to set up column numbering and communication patterns.
```
This constructor initializes the matrix structure based on the graph. The matrix values are typically set to zero. The construction process also involves determining the column numbering, distinguishing between local (diagonal block) and remote (off-diagonal block) columns.

**Distributed Assembly Workflow:**
Assembling a `DistributedCsrMatrix` involves populating its entries with values, which, like the graph construction, requires handling local and remote contributions.

1.  **`BeginAssemble()`**: This method prepares the matrix for assembly. It typically initializes all locally owned matrix entries to zero. Crucially, it also sets up communication buffers that will be used to temporarily store contributions that need to be sent to other ranks (i.e., if a local element calculation contributes to a row owned by another rank).

2.  **`Assemble()`**: Values are assembled using **global indices**.
    *   `Assemble(value, GlobalRowIndex, GlobalColIndex)`: Assembles a single scalar `value` at the specified global row and column.
    *   `Assemble(local_matrix, global_equation_ids)`: Assembles a dense `local_matrix` using a vector of `global_equation_ids`.
    *   **Local Contributions**: If `GlobalRowIndex` is owned by the current rank, the value is added to the corresponding entry in the local part of the `DistributedCsrMatrix` (either in its diagonal or off-diagonal block).
    *   **Remote Contributions**: If `GlobalRowIndex` is owned by a *different* rank, the contribution (value, global row, global column) is cached in a communication buffer. These are contributions that this rank calculated but belong to another rank's rows.
    The assembly operations themselves (the part that writes to local memory or the cache) are threadsafe.

3.  **`FinalizeAssemble()`**: This is a **collective MPI operation** and is critical.
    *   All ranks exchange the cached remote contributions. Each rank sends the contributions it has for other ranks' rows and receives contributions that other ranks have for its own rows.
    *   Received contributions are then added to the appropriate local entries.
    *   This ensures that all contributions, regardless of where they were computed, are correctly summed into the owning rank's matrix portion.
    *   All ranks must call this method concurrently.

The assembly process is illustrated below:
```cpp
    // Assuming 'A' is a DistributedCsrMatrix constructed from a dist_graph
    // and 'dist_graph' was used to define its structure.
    // Agraph in the original example is synonymous with dist_graph here.
    DistributedCsrMatrix<double, IndexType> A(dist_graph); // Construction shown above

    A.BeginAssemble();  // Prepare for assembly, init local values, setup comms

    // Loop over local elements on this rank
    for(const auto& element_info : local_elements_on_this_rank)  // connectivities in original example
    {
        LocalMatrixType K_elem; // Element's local stiffness matrix (local_matrix in original)
        EquationIdVectorType global_dofs; // Global DOF Ids for this element (equation_ids in original)

        // ... calculate K_elem and global_dofs for the current element ...

        // Assemble using global DOFs. The matrix handles local vs remote.
        A.Assemble(K_elem, global_dofs);
    }
    A.FinalizeAssemble(); // Exchange and sum remote contributions (collective call)
```

**Diagonal and Off-Diagonal Blocks:**
Each MPI rank in a `DistributedCsrMatrix` manages its locally owned rows. These rows are further conceptually divided:
*   **`DiagonalBlock`**: This is the part of the local rows where the column indices also correspond to DOFs owned by the *current* rank. It's a standard `CsrMatrix` containing purely local interactions.
*   **`OffDiagonalBlock`**: This part of the local rows contains entries where the column indices correspond to DOFs owned by *other* (remote) ranks. This is also a standard `CsrMatrix` but represents connections to ghost DOFs.

These blocks can be accessed using:
*   `GetDiagonalBlock()`: Returns a reference to a `CsrMatrix` representing the diagonal part.
*   `GetOffDiagonalBlock()`: Returns a reference to a `CsrMatrix` representing the off-diagonal part.
These are useful for certain solver algorithms or local matrix operations. The original document lists several helper functions like `GetOffDiagonalLocalIds` for mapping global to local indices within these blocks.

**Key Distributed Operations:**
*   **`SpMV(x, y)` (Sparse Matrix-Vector Product `y = A*x` or `y += A*x`)**:
    *   This is a collective operation.
    *   For the multiplication `A*x`, the vector `x` must also be distributed (`DistributedSystemVector`).
    *   Each rank computes its local part of `A*x`. This involves multiplying its `DiagonalBlock` by its local part of `x`, and its `OffDiagonalBlock` by the corresponding non-local (ghost) entries of `x`.
    *   To get these non-local entries of `x`, a `DistributedVectorImporter` object (often managed internally or provided to the `SpMV` routine) handles the necessary MPI communication to fetch ghost values from their owning ranks.
*   **`TransposeSpMV(x, y)` (Transpose SpMV `y = A^T*x` or `y += A^T*x`)**:
    *   Also a collective operation.
    *   Involves communication for contributions from the `OffDiagonalBlock` being sent back to the ranks that own those rows. A `DistributedVectorExporter` (complementary to the importer) often manages this reverse communication.

Other methods like `local_size1()` (local number of rows) and `size2()` (global number of columns) provide information about the distributed matrix dimensions. The interface for `SpMV` and `TransposeSpMV` with alpha and beta parameters is also available.

The `DistributedCsrMatrix` effectively combines the local CSR storage with the necessary communication logic (via `DataCommunicator` and internal numbering) to enable large-scale parallel matrix operations.

### Distributed System Vectors ([`DistributedSystemVector`](../Data_Structures/DistributedSystemVector.md))
A [`DistributedSystemVector`](../Data_Structures/DistributedSystemVector.md) represents a vector (like a Right-Hand Side or solution vector) partitioned across multiple MPI ranks. Each rank manages a portion of the vector's entries, corresponding to its locally owned Degrees of Freedom (DOFs).

**Role and Purpose:**
In MPI-parallel computations, global vectors (e.g., for the system `Ax=b`) must be distributed. The `DistributedSystemVector` handles this by:
*   Storing the local part of the vector on each rank.
*   Managing the `DistributedNumbering` (typically inherited from a `DistributedSparseGraph` or a `DistributedCsrMatrix`) to understand which global indices it owns.
*   Facilitating communication for assembly and operations involving non-local data access (ghost data).

**Construction:**
A `DistributedSystemVector` is typically constructed to be compatible with a `DistributedCsrMatrix`. This means it should share the same `DataCommunicator` and row partitioning (numbering).
*   **From a `DistributedSparseGraph` or `DistributedCsrMatrix`**: This is a common way to ensure compatibility. The vector will adopt the graph's or matrix's row numbering and communicator.
    ```cpp
    // Assuming 'dist_graph' is a finalized DistributedSparseGraph
    // or 'dist_matrix' is a DistributedCsrMatrix
    DistributedSystemVector<ValueType, IndexType> b(dist_graph); // or b(dist_matrix)
    // b is now partitioned consistently with the graph/matrix
    ```
*   **Explicitly with local size and communicator**: If the numbering is managed separately or needs to be defined before a matrix.
    ```cpp
    const DataCommunicator& rComm = ParallelEnvironment::GetDefaultDataCommunicator();
    IndexType local_size_for_this_rank = ...; // Derived from partitioning
    // Potentially, a full DistributedNumbering object might be passed if available
    DistributedSystemVector<ValueType, IndexType> x(local_size_for_this_rank, rComm /*, optional_numbering_info */);
    ```

**Distributed Assembly Workflow:**
Assembling contributions into a `DistributedSystemVector` (e.g., building an RHS vector `b`) mirrors the process for matrices, involving local storage and communication for remote contributions.

1.  **`BeginAssemble()`**: Prepares the vector for assembly. This typically involves:
    *   Setting all locally owned entries to zero.
    *   Initializing internal buffers for caching contributions to non-local (remote) entries.

2.  **`Assemble()` / `AssembleEntry()`**: Contributions are added using **global indices**.
    *   `Assemble(local_vector, global_equation_ids)`: Adds entries from a `local_vector` to the global vector at positions specified by `global_equation_ids`.
    *   `AssembleEntry(value, GlobalIndex)`: Adds a single `value` to the entry at `GlobalIndex`.
    *   **Local Entries**: If `GlobalIndex` is owned by the current rank, the value is added directly to the local data storage.
    *   **Non-Local Entries**: If `GlobalIndex` is owned by a *different* rank, the contribution (value, GlobalIndex) is cached.
    These assembly methods are threadsafe for concurrent calls on a single rank.

3.  **`FinalizeAssemble()`**: This is a **collective MPI operation**.
    *   The cached non-local contributions are exchanged between ranks. Kratos uses a `DistributedVectorExporter` internally for this, which efficiently sends each cached contribution to its owning rank.
    *   Each rank receives and sums up the contributions that other ranks have made to its locally owned entries.
    *   All ranks must call this method concurrently to ensure all data is correctly accumulated.

Example of assembling a distributed RHS vector `b`:
```cpp
    // Assume 'b' is a DistributedSystemVector, constructed compatibly with the system matrix
    // e.g., DistributedSystemVector<ValueType, IndexType> b(dist_graph);

    b.BeginAssemble();

    // Loop over local elements on this rank
    for(const auto& element_info : local_elements_on_this_rank)
    {
        LocalVectorType f_elem; // Element's local force vector
        EquationIdVectorType global_dofs; // Global DOF Ids for this element

        // ... calculate f_elem and global_dofs for the current element ...

        // Assemble using global DOFs. The vector handles local vs remote.
        b.Assemble(f_elem, global_dofs);
    }
    // Example of adding a single entry (e.g. a point load)
    // If specific_global_dof_id is owned by this rank, it's added locally.
    // If owned by another rank, it's cached and sent in FinalizeAssemble().
    b.AssembleEntry(load_value, specific_global_dof_id);

    b.FinalizeAssemble(); // Collective MPI call to sum all contributions
```

**Other Relevant Operations:**
Many operations are similar to the serial `SystemVector` but are executed in a distributed manner:
*   **`Norm2()` / `Norm()`**: Computes the Euclidean (L2) norm of the vector. This is a collective operation requiring MPI communication (summing partial norms from each rank).
*   **`Dot(other_distributed_vector)`**: Computes the dot product with another `DistributedSystemVector`. Also collective.
*   **`SetValue(value)`**: Sets all *locally owned* entries to a specified `value`.
*   **Local Data Access**: `operator[]` can access *locally owned* data using local indices. To work with global indices for local data, one must use the vector's `DistributedNumbering` to convert `GlobalToLocal` first. Direct access to ghost data is typically not provided through `operator[]`; ghost values are usually managed via `DistributedVectorImporter`.
*   **`CopyFrom(source_vector)`**: Copies data, potentially handling communication if `source_vector` has a different distribution or is a different type.

`DistributedSystemVector` is essential for representing global solution and RHS vectors in MPI-parallel simulations, working in tandem with `DistributedCsrMatrix` and underlying communication tools like `DistributedVectorImporter` and `DistributedVectorExporter` to manage data consistency.

### MPI Assembly Workflow Example
This section provides a conceptual C++ example demonstrating a typical distributed (MPI) workflow. It brings together the `DataCommunicator`, `DistributedSparseGraph`, `DistributedCsrMatrix`, and `DistributedSystemVector` classes to simulate a simplified Finite Element problem assembly across multiple MPI ranks.

**Note:** This code is illustrative. Actual MPI initialization and calls to get rank/size are omitted for brevity but are essential in a real MPI program. Helper functions to determine `local_size_for_this_rank` or `local_elements_on_this_rank` would depend on the specific problem and partitioning strategy.

```cpp
// Assuming necessary Kratos headers are included:
#include "includes/kratos_parameters.h"
#include "includes/parallel_environment.h"
#include "distributed_sparse_graph.h"
#include "distributed_csr_matrix.h"
#include "distributed_system_vector.h"
#include <vector>

// Define some Kratos types for clarity
using IndexType = std::size_t;
using ValueType = double;
using LocalMatrixType = Kratos::Matrix;
using LocalVectorType = Kratos::Vector;
using EquationIdVectorType = std::vector<IndexType>;
using DistributedSparseGraphType = Kratos::DistributedSparseGraph<IndexType>;
using DistributedCsrMatrixType = Kratos::DistributedCsrMatrix<ValueType, IndexType>;
using DistributedSystemVectorType = Kratos::DistributedSystemVector<ValueType, IndexType>;

void MpiAssemblyExample()
{
    // 1. Obtain DataCommunicator (represents MPI_COMM_WORLD)
    const Kratos::DataCommunicator& rComm = Kratos::ParallelEnvironment::GetDefaultDataCommunicator();
    const int current_rank = rComm.Rank();
    const int world_size = rComm.Size();

    // 2. Determine Local Problem Size (number of DOFs owned by this rank)
    // This would typically come from a distributed ModelPart after partitioning.
    // For this example, let's assume a simple division of a total number of DOFs.
    const IndexType total_global_dofs = 10; // Example total DOFs
    IndexType local_num_dofs = total_global_dofs / world_size;
    if (current_rank < total_global_dofs % world_size) {
        local_num_dofs++;
    }
    // Determine the first global DOF index owned by this rank
    IndexType first_global_dof_on_this_rank = 0;
    for(int r=0; r<current_rank; ++r) {
        IndexType dofs_on_rank_r = total_global_dofs / world_size;
        if (r < total_global_dofs % world_size) dofs_on_rank_r++;
        first_global_dof_on_this_rank += dofs_on_rank_r;
    }

    // 3. Create DistributedSparseGraph
    // Each rank specifies how many rows it owns.
    DistributedSparseGraphType dist_graph(local_num_dofs, rComm);

    // 4. Populate Graph with Connectivities (using GLOBAL indices)
    // Each rank adds entries for elements it "owns" or processes.
    // Example: Rank 0 processes an element connecting global DOFs 0, 1, (local_num_dofs_rank0 + 0)
    // Example: Rank 1 processes an element connecting global DOFs (local_num_dofs_rank0 + 0), (local_num_dofs_rank0 + 1), ...
    // For simplicity, let's assume each rank creates a few connectivities involving its own DOFs
    // and potentially DOFs from an adjacent rank (simulating a shared interface in a simple 1D decomposition).
    std::vector<EquationIdVectorType> local_element_dof_lists_global; // Store lists of global DOFs for local elements
    for(IndexType i = 0; i < local_num_dofs; ++i) {
        IndexType current_global_row = first_global_dof_on_this_rank + i;
        EquationIdVectorType dof_conn; // Connectivities for one "element" or DOF row
        dof_conn.push_back(current_global_row); // Entry for the diagonal

        // Example: connect to next global DOF if not the last one overall
        if (current_global_row + 1 < total_global_dofs) {
            dof_conn.push_back(current_global_row + 1);
        }
        // Example: connect to previous global DOF if not the first one overall
        if (current_global_row > 0) {
            dof_conn.push_back(current_global_row - 1);
        }
        local_element_dof_lists_global.push_back(dof_conn);
        dist_graph.AddEntries(dof_conn); // Add these connections to the graph
    }

    // 5. Finalize Graph (Collective MPI call)
    // Exchanges graph information between ranks to build the complete distributed structure (e.g., non-local column indices).
    dist_graph.Finalize();

    // 6. Create DistributedCsrMatrix from the graph
    DistributedCsrMatrixType A(dist_graph);

    // 7. Create DistributedSystemVectors for RHS (b) and solution (x)
    // These will be partitioned according to the graph's row distribution.
    DistributedSystemVectorType b(dist_graph);
    DistributedSystemVectorType x(dist_graph);

    // Initialize x with some values (e.g., 1.0 for all locally owned entries)
    // Note: Direct indexed access to DistributedSystemVector for writing is typically via local indices
    // or specific SetValue methods. For simplicity, we'll use SetValue for all local entries.
    x.SetValue(1.0); // Sets all *locally owned* entries on this rank to 1.0

    // 8. Assemble Global Matrix (A) and RHS vector (b)
    // This must be done on ALL ranks.
    A.BeginAssemble(); // Collective call, prepares matrix for receiving contributions
    b.BeginAssemble(); // Collective call, prepares vector for receiving contributions

    // Loop over elements processed by this rank (using the connectivities defined earlier)
    for(const auto& dof_list_global : local_element_dof_lists_global) {
        // Create dummy local stiffness matrix (LHS) and local force vector (RHS)
        // The size of local matrices/vectors depends on the number of DOFs in `dof_list_global`.
        unsigned int local_size = dof_list_global.size();
        LocalMatrixType K_local(local_size, local_size);
        LocalVectorType f_local(local_size);

        // Fill with some example values
        for(unsigned int i=0; i<local_size; ++i) {
            for(unsigned int j=0; j<local_size; ++j) {
                K_local(i,j) = (i==j) ? 2.0*local_size : -1.0;
            }
            f_local(i) = 1.0 * (dof_list_global[i] + 1); // Example value based on global DOF ID
        }

        // Assemble using GLOBAL equation IDs (dof_list_global).
        // Kratos handles whether contributions are local or need to be sent to other ranks.
        A.Assemble(K_local, dof_list_global);
        b.Assemble(f_local, dof_list_global);
    }

    // Finalize assembly (Collective MPI calls)
    // This exchanges and sums contributions that were made to non-local rows/entries.
    A.FinalizeAssemble();
    b.FinalizeAssemble();

    // 9. Perform Sparse Matrix-Vector Product (y = A*x)
    // This is also a collective operation.
    // It requires communication to handle non-local data in x (ghost values imported via DistributedVectorImporter).
    DistributedSystemVectorType y(dist_graph); // Result vector
    y.SetValue(0.0); // Ensure result vector is zero if SpMV does y += A*x
    A.SpMV(x, y);

    // Further operations (e.g., solving, output) would follow.
    // Example: Print norm of y from rank 0
    if (current_rank == 0) {
        std::cout << "MPI Workflow Example Completed. Norm of y (A*x): " << y.Norm2() << std::endl;
    }
}
```
This conceptual example outlines the primary steps in a distributed FEM assembly and matrix operation:
1.  **Initialization**: Set up MPI communication via `DataCommunicator` and determine the local portion of the problem (e.g., number of locally owned Degrees of Freedom (DOFs) and their global indices).
2.  **Graph Construction**: Each rank defines its part of the `DistributedSparseGraph` using global DOF IDs for connectivities related to its local elements/conditions. `Finalize()` is a crucial collective call to synchronize this graph structure across all ranks.
3.  **Matrix and Vector Creation**: `DistributedCsrMatrix` and `DistributedSystemVector`s are created based on the finalized distributed graph, ensuring compatible data partitioning and communication patterns.
4.  **Assembly**: Each rank computes contributions from its local elements (local stiffness matrices and force vectors) and assembles them into the global `DistributedCsrMatrix` and `DistributedSystemVector` using global DOF IDs. The `BeginAssemble()` and `FinalizeAssemble()` calls are collective and manage the necessary MPI communication to correctly sum contributions, especially those that affect DOFs (rows) owned by other ranks.
5.  **Matrix Operation**: Operations like `SpMV` (Sparse Matrix-Vector product) are performed on the distributed objects. These are also collective operations, inherently managing required data exchange (e.g., fetching non-local "ghost" values for vector `x` when computing `A*x`).

This workflow allows each MPI rank to work primarily on its subset of the problem, with Kratos's distributed linear algebra classes abstracting many of the complexities of data distribution and inter-process communication.
