---
title: SpatialSearchResultContainerVector
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: The class is designed to store and manage results from spatial searches, such as search results for multiple points.
---

# Spatial Search Result Container Vector

## Description

The class is designed to store and manage results from spatial searches, such as search results for multiple points. It provides functionalities to handle these results efficiently and supports various operations like initialization, synchronization, and result retrieval.

See [Spatial Search Result](spatial_search_result) and [Spatial Search Result Container](spatial_search_result_container).

## Implementation

The class includes error checking, ensuring that accesses are valid and that the class state is coherent before performing operations.

Can be found in [`spatial_search_result_container_vector.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/spatial_search_result_container_vector.h).

### Template arguments

- `TObjectType`: Type of the objects stored in the search results.
- `TSpatialSearchCommunication`: Defines the communication type used during spatial search (e.g., synchronous or asynchronous).

### Python exposition

Bindings live in [`add_search_strategies_to_python.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python/add_search_strategies_to_python.cpp) (`BindSpatialSearchResultContainerVector`). `SpatialSearchResultContainerVector` is mostly used internally by [`ParallelSpatialSearch`](parallel_spatial_search), which returns one to Python from every search call, but it can also be created and populated directly.

#### Instantiation

The class is exposed with four template instantiations, one per `TObjectType` (`Node`, `GeometricalObject`, `Element`, `Condition`), all of them using `SpatialSearchCommunication::SYNCHRONOUS` — `ASYNCHRONOUS` is declared in the enum but **not implemented**, and therefore not bound to Python. Each instantiation is registered under its own class name:
- `SpatialSearchResultContainerVectorNode`
- `SpatialSearchResultContainerVectorGeometricalObject`
- `SpatialSearchResultContainerVectorElement`
- `SpatialSearchResultContainerVectorCondition`

#### Methods exposed

##### Constructor (`__init__`)
Only the default (parameterless) constructor is exposed, e.g. `Kratos.SpatialSearchResultContainerVectorGeometricalObject()`.

> **Note**: `InitializeResults(NumberOfResults)` (plural, bulk pre-sizing) exists in C++ but is **not** exposed to Python — only the singular `InitializeResult()` below is. From Python, pre-size the vector by calling `InitializeResult()` once per point.

##### NumberOfSearchResults
Returns the number of search results contained in the vector.

##### InitializeResult
Initializes the search results to a default state.

##### HasResult
Checks if a result is available in the container.

##### Clear
Clears all the contents of the container, resetting it to its initial state.

##### SynchronizeAll
Synchronizes all the [`SpatialSearchResultContainer`](spatial_search_result_container) instances held by this vector across ranks in a single, batched MPI communication (building/using a shared `GlobalPointerCommunicator`), rather than synchronizing each point's container individually. This method exists **only** on `SpatialSearchResultContainerVector` — after it runs, every contained `SpatialSearchResultContainer::IsSynchronized()` becomes `true` and `GetGlobalResults()` becomes populated.

All the methods below are **bulk** queries: they exist only on `SpatialSearchResultContainerVector` (not on the individual [`SpatialSearchResultContainer`](spatial_search_result_container)), and compute the answer for every search point at once, reusing a single `GlobalPointerCommunicator` — much more efficient than looping over each point's container and resolving global pointers one by one.

##### GetDistances
Retrieves, for every search point, the list of distances to its found objects.

##### GetResultIsLocal
Retrieves, for every search point, whether each found object is local to the current rank. Takes the `DataCommunicator` as an argument.

##### GetResultRank
Retrieves, for every search point, the owning rank of each found object.

##### GetResultIsActive
Retrieves, for every search point, whether each found object is `ACTIVE`. Takes the `DataCommunicator` as an argument.

##### GetResultIsInside
Retrieves, for every search point, whether the point lies inside the geometry of each found object (only meaningful when `TObjectType` is `GeometricalObject`). Takes the points container, the `DataCommunicator`, and an optional `Tolerance` as arguments.

##### GetResultShapeFunctions
Retrieves, for every search point, the shape function values of the found object's geometry evaluated at the point (only meaningful when `TObjectType` is `GeometricalObject`). Takes the points container and the `DataCommunicator` as arguments.

##### GetResultIndices
Retrieves, for every search point, the `Id()` of each found object.

##### GetResultNodeIndices
Retrieves, for every search point, the node ids of each found object's geometry.

##### GetResultPartitionIndices
Retrieves, for every search point, the partition indices of each found object's geometry nodes.

##### GetResultCoordinates
Retrieves, for every search point, the coordinates of each found object's geometry nodes.

##### GetGlobalPointerCommunicator
Obtains the communicator used for global pointer operations across computational entities.

##### GetContainer
Accesses the underlying container (`mPointResults`) holding the search results.

#### Member variables

The class `SpatialSearchResultContainerVector` defines the following member variables:

##### `mPointResults`
This variable is a vector of pointers to `SpatialSearchResultContainerType`. It holds the results for each spatial search point.

##### `mpGlobalPointerCommunicator`
This variable is a pointer to `GlobalPointerCommunicatorType`. It facilitates communication across different data partitions, particularly useful in distributed environments.

## Example usage

### C++

#### Basic example

To help you use the `SpatialSearchResultContainerVector`, here's a simple example that demonstrates how you might use this class in your code:

```cpp
#include "includes/data_communicator.h"
#include "includes/geometrical_object.h"
#include "spatial_containers/spatial_search_result_container_vector.h"

void Example() {
    // A serial communicator works out of the box; in MPI use the model part's DataCommunicator
    DataCommunicator data_comm;

    // Create an instance of the container vector for GeometricalObjects
    SpatialSearchResultContainerVector<GeometricalObject> container_vector;

    // Pre-size the vector: one SpatialSearchResultContainer per search point
    const std::size_t number_of_points = 3;
    container_vector.InitializeResults(number_of_points);

    // Objects being "found" by the search, kept alive for the scope of the example
    std::vector<GeometricalObject> objects = {GeometricalObject(1), GeometricalObject(2), GeometricalObject(3)};

    // Add one result to each initialized position
    for (std::size_t idx = 0; idx < number_of_points; ++idx) {
        SpatialSearchResult<GeometricalObject> result(&objects[idx]);
        result.SetDistance(0.1 * idx);
        container_vector[idx].AddResult(result);
    }

    // Synchronize results across MPI ranks (a no-op with a serial communicator)
    container_vector.SynchronizeAll(data_comm);

    // Access and print out the (now synchronized) results
    for (std::size_t idx = 0; idx < number_of_points; ++idx) {
        for (auto& r_result : container_vector[idx].GetLocalResults()) {
            std::cout << "Object at index " << idx << " has distance: " << r_result.GetDistance() << std::endl;
        }
    }
}
```

### Key Points to Remember
- `InitializeResults(N)` pre-sizes the vector with `N` empty `SpatialSearchResultContainer` instances (C++ only, see the note above about the Python binding).
- `AddResult` method allows you to add spatial search results into the container at the corresponding position, obtained via `operator[]`.
- Call `SynchronizeAll` once, on the vector, after all local results have been added — this synchronizes every point in a single, batched MPI communication instead of one round-trip per point.
- Objects referenced by `SpatialSearchResult` (via `GlobalPointer`) are **not** owned by the container: their lifetime must be managed by the caller (e.g. the `ModelPart` that owns the nodes/elements/conditions), following the project's `Kratos::shared_ptr`/`unique_ptr` conventions — never `new`/`delete` them manually.

#### Advanced example

Here a code example that primarily focuses on performing spatial searches to identify and process geometrical objects relative to specified nodes within a model part. It integrates with search algorithms, data communication across computational partitions, and error-checking mechanisms.

This mirrors the pattern used by `ParallelSpatialSearch::SearchNearest` internally (see [`test_spatial_search_result_container_vector.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/tests/cpp_tests/spatial_containers/test_spatial_search_result_container_vector.cpp), `SpatialSearchResultContainerVectorGetResultShapeFunctions`, for the fully compilable version this is based on):

```cpp
// The search results, already populated by a search (e.g. SearchNearest) and containing
// one SpatialSearchResultContainer<GeometricalObject> per point in rModelPart.Nodes()
SpatialSearchResultContainerVector<GeometricalObject>& r_geometrical_object_results = ...;

// Data communicator (e.g. from the model part performing the search)
const DataCommunicator& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

// Points that were searched, in the same order used to build r_geometrical_object_results
auto& r_array_nodes = rModelPart.Nodes();

// Bulk queries: computed once for every point at once (cheaper than per-point global pointer round-trips)
const std::size_t number_of_search_results = r_geometrical_object_results.NumberOfSearchResults();
std::vector<std::vector<bool>> all_is_inside;
r_geometrical_object_results.GetResultIsInside(all_is_inside, r_array_nodes, r_data_communicator, 1.0e-5);
std::vector<std::vector<Vector>> all_shape_functions;
r_geometrical_object_results.GetResultShapeFunctions(all_shape_functions, r_array_nodes, r_data_communicator);

// Iterate over the search results, one entry per point
for (std::size_t index = 0; index < number_of_search_results; ++index) {
    auto& r_container = r_geometrical_object_results[index];

    // IsLocalPoint(): true if the point that originated this container belongs to this rank
    const bool is_local_point = r_container.IsLocalPoint();

    // At least one result found (globally) for this point
    if (r_container.NumberOfGlobalResults() > 0 && all_is_inside[index][0]) {
        const Vector& r_shape_function = all_shape_functions[index][0];

        if (is_local_point) {
            // The point belongs to this rank: link it to the found object's nodes
            for (std::size_t i = 0; i < r_shape_function.size(); ++i) {
                if (std::abs(r_shape_function[i]) > std::numeric_limits<double>::epsilon()) {
                    // Do something with the shape function (e.g. build a linear constraint)
                }
            }
        }
    }
}
```

##### How to use `r_geometrical_object_results`

1. **Retrieve search results**:
   - The number of search results is obtained using `NumberOfSearchResults()`. This number defines the loop bounds for processing each search result, and matches the number of points passed to the search.

2. **Access individual results**:
   - `r_geometrical_object_results[index]` gives the `SpatialSearchResultContainer` for the point at position `index` in the original input range.

##### Detecting local results

1. **Local vs global results**:
   - Each result container (`SpatialSearchResultContainer`) exposes `IsLocalPoint()` (does the *search point* belong to this rank?) and `NumberOfGlobalResults()`/`NumberOfLocalResults()` (how many *found objects*, globally vs. on this rank).

2. **Bulk vs per-point queries**:
   - Prefer the `SpatialSearchResultContainerVector` bulk methods (`GetResultIsInside`, `GetResultShapeFunctions`, `GetDistances`, `GetResultIndices`, ...) over calling the equivalent per-point logic manually: they compute all points in one pass over a single `GlobalPointerCommunicator`, instead of one MPI round-trip per point.

3. **Using MPI**:
   - Bulk queries already take care of the necessary MPI communication internally (through `Apply`/`GlobalPointerCommunicator`); no manual `DataCommunicator` calls are required beyond passing it in.

##### Practical steps for developers:

- **Initial setup**: Configure your environment, ensuring that the infrastructure for MPI communication is correctly set up.
- **Debugging and error checking**: Use assertions (like `KRATOS_DEBUG_ERROR_IF`) to ensure that results and operations make sense (e.g., no more than one result when expected, nodes are present in the model part).
- **Optimization**: Be mindful of performance implications, especially with MPI operations, data synchronization, and memory management (e.g., `shrink_to_fit`).

### Python

Similar to the previous C++ example, here you have a python example:

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestSpatialSearchResultContainerVector(KratosUnittest.TestCase):
    def test_vector_workflow(self):
        data_comm = Kratos.Testing.GetDefaultDataCommunicator()
        rank = data_comm.Rank()
        world_size = data_comm.Size()
        num_searches_per_rank = 2

        # 1. Create the container vector
        results_vector = Kratos.SpatialSearchResultContainerVectorGeometricalObject()

        # 2. Initialize one container per local search point
        # NOTE: only the singular InitializeResult() (no argument) is exposed to Python,
        # so pre-sizing in bulk (InitializeResults(N), C++-only) becomes a simple loop here
        for _ in range(num_searches_per_rank):
            results_vector.InitializeResult()
        self.assertEqual(results_vector.NumberOfSearchResults(), num_searches_per_rank)

        # 3. Simulate a search and add results to each container
        # NOTE: Kratos.GeometricalObject only exposes an Id-based constructor to Python
        # (unlike C++, it cannot be built directly from a Point/geometry there)
        for i in range(num_searches_per_rank):
            # Create a unique object ID based on rank and search index
            object_id = rank * num_searches_per_rank + i + 1
            found_object = Kratos.GeometricalObject(object_id)

            # Add the found object to the i-th search result container
            results_vector[i].AddResult(found_object, 0.1 * object_id)

        # 4. Synchronize ALL results across all ranks in one operation
        results_vector.SynchronizeAll(data_comm)

        # 5. Use a bulk operation to get all distances efficiently
        all_distances = results_vector.GetDistances()

        # 6. Validate the results
        total_searches = num_searches_per_rank * world_size
        self.assertEqual(len(all_distances), total_searches)

        # Each search on each rank found one object, so each sub-list has 1 distance
        for i in range(total_searches):
            self.assertEqual(len(all_distances[i]), 1)

        # Optional: Print from rank 0 for inspection
        if rank == 0:
            print(f"\nTotal number of synchronized search points: {len(all_distances)}")
            for i, dist_list in enumerate(all_distances):
                print(f"  - Search {i} found results with distances: {dist_list}")

if __name__ == '__main__':
    KratosUnittest.main()
```

#### Key Python Operations Explained

1.  **Initialization**
    * A `SpatialSearchResultContainerVector` is created, then pre-sized by calling `InitializeResult()` once per search point this rank is responsible for, appending one empty `SpatialSearchResultContainer` each time.

2.  **Populating Local Results**
    * The script loops through the initialized containers. For each one, it simulates a search by creating a sample `GeometricalObject` and adding it to the corresponding container using `results_vector[i].AddResult(...)`. At this point, all data is still local to its rank.

3.  **Bulk Synchronization**
    * A single call to `SynchronizeAll(data_comm)` performs the main communication step. It efficiently gathers all results from all sub-containers across all MPI ranks and distributes the complete global picture back to everyone.

4.  **Bulk Data Retrieval**
    * After synchronization, a method like `GetDistances()` is called. This is the key advantage of the vector class: it fetches the distances for **all** search points and **all** their results in one optimized operation, returning a list of lists.