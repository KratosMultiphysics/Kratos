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

Bindings live in [`add_search_strategies_to_python.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python/add_search_strategies_to_python.cpp) (`BindSpatialSearchResultContainer`).

#### Instantiation

The class is exposed with four template instantiations, one per `TObjectType` (`Node`, `GeometricalObject`, `Element`, `Condition`), all of them using `SpatialSearchCommunication::SYNCHRONOUS` — `ASYNCHRONOUS` is declared in the enum but **not implemented**, and therefore not bound to Python. Each instantiation is registered under its own class name:
- `SpatialSearchResultContainerNode`
- `SpatialSearchResultContainerGeometricalObject`
- `SpatialSearchResultContainerElement`
- `SpatialSearchResultContainerCondition`

#### Methods exposed

##### Constructor (`__init__`):
Only the default constructor is exposed, e.g. `Kratos.SpatialSearchResultContainerGeometricalObject()`. It does **not** take a `DataCommunicator` — synchronization is always driven from the owning [`SpatialSearchResultContainerVector`](spatial_search_result_container_vector) via `SynchronizeAll(data_communicator)`.

##### IsObjectFound
Checks if an object has been found during the spatial search.

##### NumberOfLocalResults
Returns the number of results found locally.

##### NumberOfGlobalResults
Returns the number of results found globally across all participating nodes or processes.

##### Reserve
Reserves space in the container to accommodate additional results, improving efficiency.

##### AddResult
- One version simply adds an object (`TObjectType*`) to the container.
- Another version adds an object along with its associated distance, allowing spatial contexts to consider proximity.
- When `TObjectType` is `GeometricalObject`, an additional overload accepting a full `SpatialSearchResult<GeometricalObject>` is also exposed.

##### Clear
Clears all contents of the container, resetting it to an empty state.

##### GetGlobalIndex
Retrieves the global index of an object within the container.

##### SetGlobalIndex
Sets the global index of an object within the container.

##### GetLocalIndex
Retrieves the local index of an object within the container.

##### SetLocalIndex
Sets the local index of an object within the container.

##### GetLocalResults
Retrieves a list of results that are local to the current node or process.

##### GetGlobalResults
Retrieves a list of results that are global, spanning all nodes or processes.

> **Note**: bulk-query helpers such as `GetResultIsLocal`, `GetResultRank`, `GetResultIsActive`, `GetResultIsInside`, `GetResultShapeFunctions`, `GetResultIndices`, `GetResultNodeIndices`, `GetResultPartitionIndices` and `GetResultCoordinates` are **not** methods of `SpatialSearchResultContainer` itself. They live on [`SpatialSearchResultContainerVector`](spatial_search_result_container_vector), which computes them for every point at once (more efficient than querying one `SpatialSearchResultContainer` at a time, since they share a single `GlobalPointerCommunicator`).

#### Other operations (C++ only, not exposed to Python)

##### `IsLocalPoint`
Returns `true` if `GetLocalIndex() >= 0`, i.e. if the point that originated this container's results belongs to the current rank.

##### `IsLocalSearch`
Overloaded helper that checks whether a given rank (or the rank stored at a given global index of a ranks vector) matches the current `DataCommunicator`'s rank.

##### `IsSynchronized` / `GetIsSynchronized` / `SetIsSynchronized`
Query/set whether `SynchronizeAll` has already been run for this container. `operator()`/`begin()`/`end()` (which iterate over `mGlobalResults`) raise `KRATOS_ERROR` if the container has not been synchronized yet; `operator[]` (over `mLocalResults`) works regardless.

##### `RemoveResultsFromRanksList`
Removes, from the local and global results, the entries whose rank appears in a given list of ranks. Used internally by `ParallelSpatialSearch` (`KeepOnlySmallestLambdaResult`) to discard all but the closest result once the minimum distance/rank has been determined across partitions.

#### Member variables

The class `SpatialSearchResultContainer` defines the following member variables:

##### `mLocalResults`
A vector (`std::vector`) of local search results, which stores instances of `SpatialSearchResultType` (templated by `TObjectType`). There are only the local results, until `SynchronizeAll` from `SpatialSearchResultContainerVector` is invoked.
##### `mGlobalResults`
A `GlobalPointersVector` of global search results, which stores global pointers to `SpatialSearchResultType`. This vector would be empty until `SynchronizeAll` from `SpatialSearchResultContainerVector` is invoked.
##### `mLocalIndex`
A signed index type (`std::ptrdiff_t`), used to identify local indices. Negative means that is not even a local point. This variable is added because once you pass a set of iterators to perform the search you need this information in order to identify the input point from the iterators sets from the result. This could be avoided if instead of `std::vector` we consider a map with a  certain hash, but then you need to define a hash, which could be problematic for a point which only information provided is a set of three floats with the potential issues of rounding errors. This could be avoiding using ids for nodes and hash of coordinates for the rest, but potential issues may still arise for the later.
##### `mGlobalIndex`
An unsigned index type (`std::size_t`), used for global identification purposes. At the end all the results will be globally ordered and this index helps the identifications of global results. If not defined (non existing in current rank) the value will be `0`.

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
    GeometricalObject test_object(rank + 1); // Unique ID per rank
    SpatialSearchResult<GeometricalObject> test_result(&test_object);
    test_result.SetDistance(0.5); // Set some arbitrary distance

    // Add the result to the container
    search_results.AddResult(test_result);

    // Optionally, check local results before synchronization
    auto& local_results = search_results.GetLocalResults();
    std::cout << "[Rank " << rank << "] Local results count: " << local_results.size() << std::endl;

    // NOTE: SynchronizeAll must be invoked from `SpatialSearchResultContainerVector`, so global results will be empty

    // Access synchronized results
    auto& global_results = search_results.GetGlobalResults();
    std::cout << "[Rank " << rank << "] Global results count: " << global_results.size() << std::endl; // Size will be 0 as is not synchronized

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
    - The `SpatialSearchResultContainer` is created with its (parameterless) default constructor. It does not store a `DataCommunicator` itself; synchronization across ranks is always driven externally, from the `SpatialSearchResultContainerVector` that owns it.

3. **Adding results**:
    - A `GeometricalObject` and corresponding `SpatialSearchResult` are created. The result is then added to the container. This represents a local operation.

4. **Synchronizing results**:
    - The `SynchronizeAll` method should be invoked in `SpatialSearchResultContainerVector` that holds the `SpatialSearchResultContainer` object, otherwise only local results can be accessed.

5. **Accessing results**:
    - If synchronized, you can access both local and global results through their respective methods. Otherwise only local results could be reached.

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

    # 1. Create the container
    container = Kratos.SpatialSearchResultContainerGeometricalObject()
    if rank == 0:
        print("1. Containers created on all ranks.")

    # 2. Add a unique result on each rank
    # NOTE: Kratos.GeometricalObject only exposes an Id-based constructor to Python
    # (unlike C++, it cannot be built directly from a Point/geometry there)
    object_id = rank + 1
    test_object = Kratos.GeometricalObject(object_id)
    container.AddResult(test_object)

    print(f"[Rank {rank}] Local results before sync: {container.NumberOfLocalResults()}")

    # NOTE. No global results as must be synchronized from the `SpatialSearchResultContainerVector` holder.
# To run this example, you would execute it within a Kratos MPI environment.
# run_python_example()
```

#### Key Python Operations Explained

1.  **Initialization and Setup**
    * First, we retrieve the default `DataCommunicator` and create a `ModelPart` to hold the geometric entities. The `SpatialSearchResultContainer` is then initialized with its default constructor.

2.  **Creating and Adding Local Results**
    * A unique `Node` and a corresponding `GeometricalObject` are created on each **MPI** rank. This object is then added to the container using the `AddResult` method. At this stage, each container only knows about the result created on its own rank.