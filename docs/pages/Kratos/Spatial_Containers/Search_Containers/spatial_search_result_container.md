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

Methods such as `IsObjectFound`, `Reserve`, and `Clear` provide essential utilities for result management. Additional methods like `GetDistances` and `GetResultIsInside` allow for detailed queries on the search results.

### Enumerations

- **SpatialSearchCommunication**: Defines the mode of communication for spatial search results:
  - `SYNCHRONOUS_HOMOGENEOUS`: All partitions are fully aware of the entire dataset.
  - `SYNCHRONOUS_HETEROGENEOUS`: Partitions are aware of sub-sets of the dataset.
  - `ASYNCHRONOUS`: Asynchronous communication, which is not implemented yet.

### Template arguments

- `TObjectType`: Type of the objects stored in the search results.
- `TSpatialSearchCommunication`: Defines the communication type used during spatial search (e.g., synchronous or asynchronous). See previous point.

### Python exposition

#### Instantiation

The class is exposed with several template instantiations for different types (`Node`, `GeometricalObject`, `Element`, `Condition`) with two communication types:
- **SYNCHRONOUS_HOMOGENEOUS**: This mode presumably ensures a consistent, uniform handling of communication across all nodes or elements during spatial searches, regardless of their heterogeneity or distribution.
- **SYNCHRONOUS_HETEROGENEOUS**: In contrast, this mode might allow for different handling or processing strategies depending on the node or element type or state, providing more flexibility in operations that require acknowledgment of diversity in the dataset or environment.

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

##### GetDistances
Returns the list of distances associated with the objects in the container.

See [Spatial Search Result Container Vector GetDistances](spatial_search_result_container_vector#GetDistances).

##### GetResultIsLocal
Checks if a result is local to the current process or node.

See [Spatial Search Result Container Vector GetResultIsLocal](spatial_search_result_container_vector#GetResultIsLocal).

##### GetResultRank
Retrieves the rank of the process or node where a result was found.

See [Spatial Search Result Container Vector GetResultRank](spatial_search_result_container_vector#GetResultRank).

##### GetResultIsActive
Checks if the result is considered active based on the spatial search criteria.

See [Spatial Search Result Container Vector GetResultIsActive](spatial_search_result_container_vector#GetResultIsActive).

##### GetResultIsInside
Determines whether the result is inside a specified boundary or region.

See [Spatial Search Result Container Vector GetResultIsInside](spatial_search_result_container_vector#GetResultIsInside).

##### GetResultShapeFunctions
Retrieves the shape functions associated with the results, useful in finite element methods.

See [Spatial Search Result Container Vector GetResultShapeFunctions](spatial_search_result_container_vector#GetResultShapeFunctions).

##### GetResultIndices
Retrieves indices associated with the results, useful for referencing within larger datasets or arrays.

See [Spatial Search Result Container Vector GetResultIndices](spatial_search_result_container_vector#GetResultIndices).

##### GetResultNodeIndices
Retrieves the node indices specifically, used in mesh-based applications.

See [Spatial Search Result Container Vector GetResultNodeIndices](spatial_search_result_container_vector#GetResultNodeIndices).

##### GetResultPartitionIndices
Retrieves indices useful for understanding data partitioning in distributed environments.

See [Spatial Search Result Container Vector GetResultPartitionIndices](spatial_search_result_container_vector#GetResultPartitionIndices).

##### GetResultCoordinates
Gets the coordinates of the result, useful in spatial applications.

See [Spatial Search Result Container Vector GetResultCoordinates](spatial_search_result_container_vector#GetResultCoordinates).

##### GetDataCommunicator
Accesses the `DataCommunicator` used for data synchronization and communication.

See [Spatial Search Result Container Vector GetDataCommunicator](spatial_search_result_container_vector#GetDataCommunicator).

##### GetLocalResults
Retrieves a list of results that are local to the current node or process.

##### GetGlobalResults
Retrieves a list of results that are global, spanning all nodes or processes.

##### GetGlobalPointerCommunicator
Provides access to the global pointer communicator, essential for managing pointers in a distributed system.

#### Member variables

The class `SpatialSearchResultContainer` defines the following member variables:

##### `mrDataCommunicator`
A reference to the `DataCommunicator` used for communication between different data processors.
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
#include "spatial_containers/spatial_search_result_container.h"

int main() {
    // Create a data communicator
    DataCommunicator& r_data_comm = DataCommunicator::GetDefault();

    // Create a container for storing search results
    SpatialSearchResultContainer<GeometricalObject> search_results(r_data_comm);

    // Create a geometric object and a search result
    GeometricalObject test_object = GeometricalObject(r_data_comm.Rank() + 1); // Unique ID per rank
    SpatialSearchResult<GeometricalObject> test_result(&test_object);
    test_result.SetDistance(0.5); // Set some arbitrary distance

    // Add the result to the container
    search_results.AddResult(test_result);

    // Optionally, check local results before synchronization
    auto& local_results = search_results.GetLocalResults();
    std::cout << "Local results count: " << local_results.size() << std::endl;

    // Synchronize results across all MPI nodes
    search_results.SynchronizeAll();

    // Access synchronized results
    auto& global_results = search_results.GetGlobalResults();
    std::cout << "Global results count: " << global_results.size() << std::endl;

    // Clear the results in the container
    search_results.Clear();
    std::cout << "Results cleared. Local results count: " << search_results.NumberOfLocalResults() << std::endl;

    return 0;
}
```

### Key operations explained:

1. **Initialization and data communication**:
    - Start by initializing MPI and retrieving the default data communicator, which handles data transfer between different computational nodes.

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