---
title: SystemVector
keywords: SystemVector
tags: [SystemVector, Data Structures, Linear Algebra]
sidebar: kratos_for_developers
summary: Comprehensive documentation for the Kratos SystemVector class, system of equations vector class.
---

## SystemVector Class

The `SystemVector` is a sequential data structure used in the Finite Element Method (FEM) for assembling and managing system vectors. These vectors typically store solution unknowns or right-hand side contributions to the global system of equations. It provides a convenient interface for common vector operations and assembly processes within a serial (non-distributed) context.

### Constructors

The `SystemVector` class provides several constructors for flexibility:

*   **Default Constructor:**
    ```cpp
    SystemVector() = default;
    ```
    Creates an empty `SystemVector`.

*   **From Graph (Size Initialization):**
    ```cpp
    SystemVector(const SparseGraph<IndexType>& rGraph);
    SystemVector(const SparseContiguousRowGraph<IndexType>& rGraph);
    ```
    Constructs a `SystemVector` with a size determined by the provided `SparseGraph` or `SparseContiguousRowGraph`. The underlying data communicator is also taken from the graph.

*   **With Size and Communicator:**
    ```cpp
    SystemVector(IndexType size, DataCommunicator& rComm = ParallelEnvironment::GetDataCommunicator("Serial"));
    ```
    Creates a `SystemVector` of a specified `size`. An optional `DataCommunicator` can be provided; if not, it defaults to a serial communicator. This constructor will raise an error if a distributed communicator is passed, as `SystemVector` is a serial (non-distributed) container.

*   **From `Vector` and Communicator:**
    ```cpp
    SystemVector(const Vector& data, DataCommunicator& rComm = ParallelEnvironment::GetDataCommunicator("Serial"));
    ```
    Constructs a `SystemVector` by copying data from an existing Kratos `Vector`. Similar to the size-based constructor, it accepts an optional `DataCommunicator` and defaults to serial. It will also error if a distributed communicator is provided.

*   **Copy Constructor:**
    ```cpp
    explicit SystemVector(const SystemVector<TDataType,TIndexType>& rOtherVector);
    ```
    Creates a deep copy of another `SystemVector`.

*   **Move Constructor:**
    ```cpp
    SystemVector(SystemVector<TDataType,TIndexType>&& rOtherVector);
    ```
    Moves the resources from another `SystemVector`, leaving the other vector in a valid but unspecified state.

### Operators

`SystemVector` supports common arithmetic and access operators:

*   **Assignment (`=`):**
    ```cpp
    SystemVector& operator=(SystemVector const& rOtherVector); // Copy assignment
    SystemVector& operator=(SystemVector&& rOtherVector);      // Move assignment
    ```
    Assigns the content of another `SystemVector`.

*   **Addition Assignment (`+=`):**
    ```cpp
    SystemVector& operator+=(const SystemVector& rOtherVector);
    ```
    Adds another `SystemVector` element-wise.

*   **Subtraction Assignment (`-=`):**
    ```cpp
    SystemVector& operator-=(const SystemVector& rOtherVector);
    ```
    Subtracts another `SystemVector` element-wise.

*   **Scalar Multiplication Assignment (`*=`):**
    ```cpp
    SystemVector& operator*=(const TDataType multiplier_factor);
    ```
    Multiplies each element by a scalar.

*   **Scalar Division Assignment (`/=`):**
    ```cpp
    SystemVector& operator/=(const TDataType divide_factor);
    ```
    Divides each element by a scalar.

*   **Element Access (`()` and `[]`):**
    ```cpp
    TDataType& operator()(IndexType I);
    const TDataType& operator()(IndexType I) const;
    TDataType& operator[](IndexType I);
    const TDataType& operator[](IndexType I) const;
    ```
    Provides read and write access to individual elements using their index.

### Key Methods

Here are some of the main methods provided by `SystemVector`:

*   **`Clear()`**
    ```cpp
    void Clear();
    ```
    Removes all elements from the vector, resulting in a size of 0.

*   **`SetValue(const TDataType value)`**
    ```cpp
    void SetValue(const TDataType value);
    ```
    Assigns the given `value` to all elements in the vector.

*   **`size() const`**
    ```cpp
    IndexType size() const;
    ```
    Returns the number of elements currently stored in the vector.

*   **`data()`**
    ```cpp
    DenseVector<TDataType>& data();
    const DenseVector<TDataType>& data() const;
    ```
    Returns a reference (or const reference) to the underlying `DenseVector<TDataType>` that stores the vector's data. This allows for low-level access if needed.

*   **`Add(const TDataType factor, const SystemVector& rOtherVector)`**
    ```cpp
    void Add(const TDataType factor, const SystemVector& rOtherVector);
    ```
    Performs the operation `this = this + factor * rOtherVector`.

*   **`Norm() const`**
    ```cpp
    TDataType Norm() const;
    ```
    Calculates and returns the L2 (Euclidean) norm of the vector.

*   **`Dot(const SystemVector& rOtherVector, IndexType gather_on_rank=0)`**
    ```cpp
    TDataType Dot(const SystemVector& rOtherVector, IndexType gather_on_rank=0);
    ```
    Computes the dot product between this vector and `rOtherVector`. The `gather_on_rank` parameter is included for interface consistency with `DistributedSystemVector` and is not used in `SystemVector`.

*   **`BeginAssemble()` and `FinalizeAssemble()`**
    ```cpp
    void BeginAssemble();
    void FinalizeAssemble();
    ```
    These methods are part of the assembly interface. In the serial `SystemVector`, they currently perform no operation. They are included to maintain a consistent interface with the `DistributedSystemVector` where they manage communication during parallel assembly.

*   **`Assemble(const TVectorType& rVectorInput, const TIndexVectorType& EquationId)`**
    ```cpp
    template<class TVectorType, class TIndexVectorType >
    void Assemble(
        const TVectorType& rVectorInput,
        const TIndexVectorType& EquationId
    );
    ```
    Assembles contributions from a local `rVectorInput` into the `SystemVector` at the global indices specified by `EquationId`. This is a key method for FEM assembly. `AtomicAdd` is used internally for thread-safe assembly.

*   **`AssembleEntry(const TDataType rValue, const IndexType GlobalI)`**
    ```cpp
    void AssembleEntry(
        const TDataType rValue,
        const IndexType GlobalI
    );
    ```
    Atomically adds `rValue` to the element at the specified `GlobalI` in the vector.

*   **`GetComm() const`**
    ```cpp
    const DataCommunicator& GetComm() const;
    ```
    Returns a const reference to the `DataCommunicator` associated with the vector. For `SystemVector`, this is typically the serial communicator.

*   **`PrintInfo(std::ostream& rOStream) const` and `PrintData(std::ostream& rOStream) const`**
    ```cpp
    void PrintInfo(std::ostream& rOStream) const;
    void PrintData(std::ostream& rOStream) const;
    ```
    `PrintInfo` outputs a general description of the `SystemVector` object, while `PrintData` prints the actual data contained within the vector to the provided output stream.

### Code Examples

Below are examples demonstrating the usage of `SystemVector`.

#### **1. Construction and Initialization**

```cpp
#include "containers/system_vector.h"
#include "includes/data_communicator.h"
#include "includes/parallel_environment.h"

// ...

// Default constructor
Kratos::SystemVector<> vec1;

// Constructor with size (using serial communicator by default)
Kratos::SystemVector<> vec2(10); // Creates a vector of size 10

// Initialize all elements to a value
vec2.SetValue(5.0);

// Constructor from a Kratos Vector
Kratos::Vector kratos_vec(5, 1.23); // A standard Kratos Vector
Kratos::SystemVector<> vec3(kratos_vec);

// Copy constructor
Kratos::SystemVector<> vec4(vec2);

// Using a graph to define size
// Assuming 'graph' is a SparseContiguousRowGraph or SparseGraph object
Kratos::SystemVector<> vec_from_graph(graph);
```
#### **2. Common Operations**

```cpp
#include "containers/system_vector.h"
#include "includes/data_communicator.h" // For DataCommunicator
#include "includes/parallel_environment.h" // For ParallelEnvironment
// ...
Kratos::SystemVector<double> a(4);
a.SetValue(5.0); // a = [5.0, 5.0, 5.0, 5.0]
Kratos::SystemVector<double> b(4);
b.SetValue(3.0); // b = [3.0, 3.0, 3.0, 3.0]
// Element access (0-based indexing)
a[0] = 10.0;       // a = [10.0, 5.0, 5.0, 5.0]
double val = b(1); // val = 3.0
// Arithmetic operations
Kratos::SystemVector<double> c = a; // c is a copy of a
c += b;                             // c = [13.0, 8.0, 8.0, 8.0]
c -= b;                             // c = [10.0, 5.0, 5.0, 5.0] (back to a's state after modification)
c.Add(2.0, b);                      // c = c + 2.0 * b = [10+6, 5+6, 5+6, 5+6] = [16.0, 11.0, 11.0, 11.0]
c *= 0.5;                           // c = [8.0, 5.5, 5.5, 5.5]
c /= 2.0;                           // c = [4.0, 2.75, 2.75, 2.75]
// Norm
double norm_a = a.Norm(); // Calculates L2 norm of a
// Dot product
double dot_ab = a.Dot(b);
// Get underlying data (advanced)
Kratos::DenseVector<double>& raw_data_a = a.data();
raw_data_a[1] = 7.0; // Modifies a directly: a = [10.0, 7.0, 5.0, 5.0]
// Print
a.PrintInfo(std::cout);
a.PrintData(std::cout);
```

#### **3. Assembly Process**

This example mimics the assembly typically done in an FEM context.

```cpp
#include "containers/system_vector.h"
#include "containers/sparse_contiguous_row_graph.h" // For graph example
#include "includes/element.h" // For connectivities (conceptual)

// ...

// Assume 'connectivities' is a std::vector of std::vector<int>,
// representing element node indices (Equation IDs).
// Example: Element 0 connects DOFs {0, 1}, Element 1 connects DOFs {1, 2, 3}
std::vector<std::vector<Kratos::SystemVector<>::IndexType>> connectivities;
connectivities.push_back({0, 1});
connectivities.push_back({1, 2, 3});
connectivities.push_back({0, 3, 4});
// ... add more element connectivities

// Create a graph to determine the size of the system vector
// In a real scenario, this graph would be built based on the mesh.
Kratos::SparseContiguousRowGraph<> fem_graph(5); // Assuming 5 DOFs in total for this example
// Populate the graph (omitted for brevity, but AddEntries would be used)
// For SystemVector, the graph primarily defines the size.

// Create the system vector (e.g., for the RHS or solution vector)
Kratos::SystemVector<> b_rhs(fem_graph);
b_rhs.SetValue(0.0); // Initialize to zero

// Simulate assembling contributions from elements
b_rhs.BeginAssemble(); // No-op in SystemVector, but good practice

// Contribution from Element 0
Kratos::Vector element0_contribution(2);
element0_contribution[0] = 1.5; // Contribution to DOF 0
element0_contribution[1] = 0.5; // Contribution to DOF 1
b_rhs.Assemble(element0_contribution, connectivities[0]);
// b_rhs is now approximately: [1.5, 0.5, 0.0, 0.0, 0.0] (depending on graph finalization)

// Contribution from Element 1
Kratos::Vector element1_contribution(3);
element1_contribution[0] = 0.2; // Contribution to DOF 1
element1_contribution[1] = 0.8; // Contribution to DOF 2
element1_contribution[2] = 0.3; // Contribution to DOF 3
b_rhs.Assemble(element1_contribution, connectivities[1]);
// b_rhs is now approximately: [1.5, 0.5+0.2, 0.8, 0.3, 0.0] = [1.5, 0.7, 0.8, 0.3, 0.0]

// Assemble a single entry
b_rhs.AssembleEntry(10.0, 4); // Add 10.0 to DOF 4
// b_rhs is now approximately: [1.5, 0.7, 0.8, 0.3, 10.0]


b_rhs.FinalizeAssemble(); // No-op in SystemVector

// b_rhs now contains the assembled values.
b_rhs.PrintData(std::cout);
```