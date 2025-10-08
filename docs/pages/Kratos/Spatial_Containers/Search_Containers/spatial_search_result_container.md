---
title: SpatialSearchResultContainer
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: A templated class designed to store and manage the results of spatial searches.
---

# Spatial Search Result Container

## Description

This is a templated templated class designed to store and manage the results of spatial searches. Specially designed for the management and communication of spatial search results across different computational partitions in a distributed computing environment.

See [Spatial Search Result](spatial_search_result) and [Spatial Search Result Container Vector](spatial_search_result_container_vector).

## Implementation

Can be found in [`spatial_search_result_container.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/spatial_search_result_container.h).

Utilizes global pointers to manage data references across different computational nodes, enhancing the robustness and flexibility of distributed data handling.

Methods such as `IsObjectFound`, `Reserve`, and `Clear` provide essential utilities for result management.

### Enumerations

- **`SpatialSearchCommunication`**: Defines the communication mode for spatial search results:
    - `SYNCHRONOUS`: All partitions gather results from all other partitions to build a complete global list.
    - `ASYNCHRONOUS`: Asynchronous communication, which is **not implemented yet**.

### Template arguments

- `TObjectType`: Type of the objects stored in the search results.
- `TSpatialSearchCommunication`: Defines the communication type used during spatial search (e.g., synchronous or asynchronous). See previous point.

### Python exposition

#### Instantiation

The class is exposed with several template instantiations for different types (`Node`, `GeometricalObject`, `Element`, `Condition`) with two communication types:
- **SYNCHRONOUS**: This mode presumably ensures a consistent, uniform handling of communication across all nodes or elements during spatial searches, regardless of their heterogeneity or distribution.
- **ASYNCHRONOUS**: Asynchronous communication, which is **not implemented yet**.

Each type and communication mode combination is then bound to the module with a specific class name, such as `SpatialSearchResultContainerNode` for nodes with homogeneous communication, or `SpatialSearchResultContainerNodeHeterogeneous` for nodes with heterogeneous communication.

#### Methods exposed

##### Constructor (`__init__`):
Initializes the container with a reference to a `DataCommunicator` instance.

##### IsObjectFound
Checks if an object has been found during the spatial search.

##### NumberOfLocalResults
Returns the number of results found locally.

##### NumberOfGlobalResults
Returns the number of results found globally across all participating nodes or processes.

##### Reserve
Reserves space in the container to accommodate additional results, improving efficiency.

##### AddResult
- One version simply adds an object to the container.
- Another version adds an object along with its associated distance, allowing spatial contexts to consider proximity.

##### Clear
Clears all contents of the container, resetting it to an empty state.

##### SynchronizeAll
Synchronizes the container's contents across all nodes or processes involved in the spatial search, ensuring consistency.

See [Spatial Search Result Container Vector SynchronizeAll](spatial_search_result_container_vector#SynchronizeAll).

##### GetGlobalIndex
Retrieves the global index of an object within the container.

##### SetGlobalIndex
Sets the global index of an object within the container.

##### GetLocalIndex
Retrieves the local index of an object within the container.

##### SetLocalIndex
Sets the local index of an object within the container.

##### GetResultIsLocal
Checks if a result is local to the current process or node.

See [Spatial Search Result Container Vector GetResultIsLocal](spatial_search_result_container_vector#GetResultIsLocal).

##### GetLocalResults
Retrieves a list of results that are local to the current node or process.

##### GetGlobalResults
Retrieves a list of results that are global, spanning all nodes or processes.

##### GetGlobalPointerCommunicator
Provides access to the global pointer communicator, essential for managing pointers in a distributed system.

#### Member variables

The class `SpatialSearchResultContainer` defines the following member variables:

##### `mLocalResults`
A vector (`std::vector`) of local search results, which stores instances of `SpatialSearchResultType` (templated by `TObjectType`). There are only the local results, until `SynchronizeAll` is invoked.
##### `mGlobalResults`
A `GlobalPointersVector` of global search results, which stores global pointers to `SpatialSearchResultType`. This vector would be empty until `SynchronizeAll` is invoked.
##### `mLocalIndex`
A signed index type (`std::ptrdiff_t`), used to identify local indices. Negative means that is not even a local point. This variable is added because once you pass a set of iterators to perform the search you need this information in order to identify the input point from the iterators sets from the result. This could be avoided if instead of `std::vector` we consider a map with a  certain hash, but then you need to define a hash, which could be problematic for a point which only information provided is a set of three floats with the potential issues of rounding errors. This could be avoiding using ids for nodes and hash of coordinates for the rest, but potential issues may still arise for the later.
##### `mGlobalIndex`
An unsigned index type (`std::size_t`), used for global identification purposes. At the end all the results will be globally ordered and this index helps the identifications of global results. If not defined (non existing in current rank) the value will be `0`.
##### `mpGlobalPointerCommunicator`
A pointer to a `GlobalPointerCommunicator` of `SpatialSearchResultType`, used to manage global data communication.

## Example usage

### C++

To help you use the `SpatialSearchResultContainer`, here's a simple example that demonstrates how you might use this class in your code:

```cpp
#include <iostream>
#include "includes/data_communicator.h"
#include "includes/node.h"
#include "spatial_containers/spatial_search_result_container.h"

int main() {
    // Create a data communicator
    const DataCommunicator& r_data_comm = DataCommunicator::GetDefault();
    const int rank = r_data_comm.Rank();

    // Create a container for storing search results using the default constructor
    SpatialSearchResultContainer<GeometricalObject> search_results;

    // Create a geometric object and a search result in a more complete way
    Node::Pointer p_node = Kratos::make_intrusive<Node>(rank + 1, 0.0, 0.0, 0.0); // Unique ID per rank
    GeometricalObject test_object(Kratos::Point(p_node));
    SpatialSearchResult<GeometricalObject> test_result(&test_object);
    test_result.SetDistance(0.5); // Set some arbitrary distance

    // Add the result to the container
    search_results.AddResult(test_result);

    // Optionally, check local results before synchronization
    auto& local_results = search_results.GetLocalResults();
    std::cout << "[Rank " << rank << "] Local results count: " << local_results.size() << std::endl;

    // Synchronize results across all MPI nodes, passing the communicator
    search_results.SynchronizeAll(r_data_comm);

    // Access synchronized results
    auto& global_results = search_results.GetGlobalResults();
    std::cout << "[Rank " << rank << "] Global results count: " << global_results.size() << std::endl;

    // Clear the results in the container
    search_results.Clear();
    std::cout << "[Rank " << rank << "] Results cleared. Local results count: " << search_results.NumberOfLocalResults() << std::endl;

    return 0;
}
```

### Key operations explained:

1. **Initialization and data communication**:
    - Start by initializing **MPI** and retrieving the default data communicator, which handles data transfer between different computational nodes.

2. **Creating the container**:
    - The `SpatialSearchResultContainer` is initialized using the data communicator to handle the synchronization of data across different processes.

3. **Adding results**:
    - A `GeometricalObject` and corresponding `SpatialSearchResult` are created. The result is then added to the container. This represents a local operation.

4. **Synchronizing results**:
    - The `SynchronizeAll` method is called to sync results across all processes. This ensures that all nodes have the same global view of the search results.

5. **Accessing results**:
    - After synchronization, you can access both local and global results through their respective methods.

6. **Clearing results**:
    - The `Clear` method resets the container, removing all stored results, both locally and globally, if synchronized.

This example code illustrates a typical usage scenario for managing spatial search results in a distributed system using the `SpatialSearchResultContainer`.

### Python

Similar to the previous C++ example, here you have a python example:

```python
import KratosMultiphysics as Kratos

def run_python_example():
    """Demonstrates the SpatialSearchResultContainer workflow in Python."""
    data_comm = Kratos.Testing.GetDefaultDataCommunicator()
    rank = data_comm.Rank()
    world_size = data_comm.Size()

    # In a real case, the ModelPart would be part of a larger analysis
    model = Kratos.Model()
    model_part = model.CreateModelPart("Main")

    # 1. Create the container
    container = Kratos.SpatialSearchResultContainerGeometricalObject()
    if rank == 0:
        print("1. Containers created on all ranks.")

    # 2. Add a unique result on each rank
    node_id = rank + 1
    new_node = model_part.CreateNewNode(node_id, 0.0, 0.0, 0.0)
    test_object = Kratos.GeometricalObject(Kratos.Point2D(new_node))
    container.AddResult(test_object)

    print(f"[Rank {rank}] Local results before sync: {container.NumberOfLocalResults()}")

    # 3. Synchronize results across all ranks
    # Note: No explicit barrier is needed; the collective call handles synchronization.
    if rank == 0:
        print("\n2. Synchronizing all containers...")
    container.SynchronizeAll(data_comm)

    # 4. Access the global results (output from rank 0 for clarity)
    if rank == 0:
        print("\n3. Synchronization Complete.")
        num_global = container.NumberOfGlobalResults()
        print(f"   Total global results: {num_global}")

        # In a test, we would assert this
        assert num_global == world_size

        for result_ptr in container: # Iterating over global results
            print(f"     - Found result from Rank {result_ptr.Rank()}" \
                  f" (Object ID: {result_ptr.Get().Id})")

# To run this example, you would execute it within a Kratos MPI environment.
# run_python_example()
```

#### Key Python Operations Explained

1.  **Initialization and Setup**
    * First, we retrieve the default `DataCommunicator` and create a `ModelPart` to hold the geometric entities. The `SpatialSearchResultContainer` is then initialized with its default constructor.

2.  **Creating and Adding Local Results**
    * A unique `Node` and a corresponding `GeometricalObject` are created on each **MPI** rank. This object is then added to the container using the `AddResult` method. At this stage, each container only knows about the result created on its own rank.

3.  **Synchronization**
    * The `SynchronizeAll` method is called, passing the `DataCommunicator`. This is the key step where the container communicates with all other ranks, gathering all the individual local results into one globally consistent list.

4.  **Accessing Global Results**
    * After synchronization, the container can be iterated over directly (e.g., `for result in container:`). This provides access to the complete list of results from **all** ranks, enabling global operations or analysis.

5.  **Clearing the Container**
    * The `Clear` method can be called to empty the container of both its local and global results, preparing it for a new search operation.