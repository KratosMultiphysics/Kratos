---
title: ParallelSpatialSearch
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: This is a parallel spatial search ready for MPI searches
---

# Parallel Spatial Search

## Description

`ParallelSpatialSearch` is a generic search wrapper that must be adapted and specialized for every search object (backend). It is designed to efficiently handle searches across different processes in a distributed computing environment using MPI, but it also works transparently in serial runs (with a serial `DataCommunicator`).

Internally it takes any of the *Kratos* spatial containers (see [`GeometricalObjectsBins`](geometrical_objects_bins), [KDTree](kd_tree), [OCTree](octree), [`BinsStatic`](bins_static) and [`BinsDynamic`](bins_dynamic)) as its local (per-rank) search backend, and adds the logic required to synchronize points and results across partitions: it synchronizes the points to be searched (together with the local bounding boxes) to every rank that might contain a result, performs the local search on each rank in parallel (via `IndexPartition`), and finally synchronizes the results back with [`SpatialSearchResultContainerVector::SynchronizeAll`](spatial_search_result_container_vector#SynchronizeAll).

See [`SpatialSearchResultContainer`](spatial_search_result_container) and [`SpatialSearchResultContainerVector`](spatial_search_result_container_vector), which are used to store and synchronize the results.

## Implementation

Can be found in [`parallel_spatial_search.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/parallel_spatial_search.h) / [`parallel_spatial_search.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/parallel_spatial_search.cpp).

### Template arguments

- `TSearchObject`: The type of the local search backend, e.g. `GeometricalObjectsBins` or one of the `Tree<...>`/`BinsDynamic<...>` instantiations. It must expose `ObjectType`, `PointType` and a constructor compatible with the backend (bins/tree/`bucket_size`).
- `TSpatialSearchCommunication`: The communication strategy. It considers the enum `SpatialSearchCommunication`, see [`SpatialSearchResultContainer`](spatial_search_result_container):
  - `SYNCHRONOUS`: All partitions gather the results from all other partitions, so every rank ends up with a complete, globally synchronized view of the results. This is the only strategy currently implemented, and the only one exposed to Python.
  - `ASYNCHRONOUS`: Asynchronous communication, **which is not implemented yet**.

### Construction

The object to be searched (nodes, elements or conditions) is passed once at construction time and cannot be changed afterwards — a new `ParallelSpatialSearch` instance must be created if the underlying entities change. Internally, for `GeometricalObjectsBins` the geometrical objects are stored directly in the bins; for the tree/bins-based backends, the entities are first converted to a `PointVector` (`SearchUtilities::PreparePointsSearch`) before building the search structure. `Parameters` (see [Parameters & Defaults](#parameters--defaults)) can optionally be passed to configure `allocation_size` and `bucket_size`.

### Search operations

All search operations take a range of point iterators (e.g. `ModelPart::NodesContainerType` iterators) and populate a [`SpatialSearchResultContainerVector`](spatial_search_result_container_vector), one entry per input point, in the same order as the input range:

- `SearchInRadius`: Finds every object within a given radius of each point. Results include the distance to each object found.
- `SearchNearestInRadius`: Finds the single nearest object within a given radius of each point. If several objects share the minimum distance only one is kept; if none is found the result is marked as not found.
- `SearchNearest`: Like `SearchNearestInRadius`, but with no radius limit (the global bounding box diagonal is used internally as the search radius).
- `SearchIsInside`: Checks whether each point lies inside one of the geometrical objects of the domain (only available when `TSearchObject` is `GeometricalObjectsBins`, see `IsGeometricalObjectBins`).

`GetGlobalBoundingBox` returns the bounding box that encloses the local bounding boxes of all the participating ranks.

## Python exposition

Bindings live in [`add_search_strategies_to_python.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python/add_search_strategies_to_python.cpp). Only the `SYNCHRONOUS` communication type is currently exposed to Python.

#### Exposed methods

1. **Constructors**:
   - Different constructors are exposed depending on `ObjectType` (`TSearchObject::ObjectType`):
     - `Node` → `(NodesContainerType&, DataCommunicator&)` and `(NodesContainerType&, DataCommunicator&, Parameters)`.
     - `Element` (or `GeometricalObject`) → the equivalent constructors taking `ElementsContainerType&`.
     - `Condition` (or `GeometricalObject`) → the equivalent constructors taking `ConditionsContainerType&`.
   - Since `GeometricalObjectsBins::ObjectType` is `GeometricalObject`, `ParallelSpatialSearchGeometricalObjectBins` exposes both the `Elements` and `Conditions` constructors.

2. **Spatial Search Methods** (all take a `NodesContainerType`, i.e. `model_part.Nodes` or `sub_model_part.Nodes`, and return a [`SpatialSearchResultContainerVector`](spatial_search_result_container_vector)):
   - `GetGlobalBoundingBox()`: Retrieves the global bounding box of the search space.
   - `SearchInRadius(nodes, radius)`: Performs a radius-based search.
   - `SearchNearestInRadius(nodes, radius)`: Searches for the nearest object within a specified radius.
   - `SearchNearest(nodes)`: Searches for the nearest object.
   - `SearchIsInside(nodes)`: Checks if the points are inside a geometrical object (only exposed when `TSearchObject` is `GeometricalObjectsBins`).

#### Instantiation

Each combination of backend and object type is bound under its own Python class name:

##### 1. `ParallelSpatialSearchGeometricalObjectBins`
   - Wraps [`GeometricalObjectsBins`](geometrical_objects_bins).
   - Exposes constructors for `Elements` and `Conditions`, plus all four spatial search methods (including `SearchIsInside`).

##### 2. KDTree wrappers
   - `ParallelSpatialSearchKDTreeNode`, `ParallelSpatialSearchKDTreeElement`, `ParallelSpatialSearchKDTreeCondition`.

##### 3. OCTree wrappers
   - `ParallelSpatialSearchOCTreeNode`, `ParallelSpatialSearchOCTreeElement`, `ParallelSpatialSearchOCTreeCondition`.

##### 4. Static and dynamic bins wrappers
   - Static: `ParallelSpatialSearchStaticBinsTreeNode`, `ParallelSpatialSearchStaticBinsTreeElement`, `ParallelSpatialSearchStaticBinsTreeCondition`.
   - Dynamic: `ParallelSpatialSearchDynamicBinsNode`, `ParallelSpatialSearchDynamicBinsElement`, `ParallelSpatialSearchDynamicBinsCondition`.

The KDTree, OCTree and static bins wrappers do **not** expose `SearchIsInside`, since that operation is only meaningful (and implemented) for `GeometricalObjectsBins`.

## Example usage

### C++

Here's an example of how you can utilize `ParallelSpatialSearch<GeometricalObjectsBins>` in C++ (see [`test_parallel_spatial_search.cpp`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/tests/cpp_tests/spatial_containers/test_parallel_spatial_search.cpp) for more examples, including tree/bins backends):

```cpp
#include "containers/model.h"
#include "tests/test_utilities/cpp_tests_utilities.h"
#include "spatial_containers/parallel_spatial_search.h"

using namespace Kratos;

// A convenience alias, as done in the Kratos test suite
using ParallelSpatialSearchGeometricalObjectsBins = ParallelSpatialSearch<GeometricalObjectsBins>;

static void Example() {
    Model current_model;

    // Create a cube skin model part (the objects to be searched)
    ModelPart& cube_skin = CppTestsUtilities::CreateCubeSkinModelPart(current_model);

    // A serial data communicator works out of the box; in MPI use the model part's DataCommunicator
    DataCommunicator serial_communicator;

    // Instantiate the parallel spatial search with the elements of the model part
    ParallelSpatialSearchGeometricalObjectsBins parallel_spatial_search(cube_skin.Elements(), serial_communicator);

    // Create a model part for the points to be searched
    ModelPart& point_model_part = current_model.CreateModelPart("PointModelPart");
    point_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    // Create a vector to store the search results (one entry per input point)
    ParallelSpatialSearchGeometricalObjectsBins::ResultContainerVectorType results;

    // Conduct a search in radius of 0.3 units around the newly added node
    parallel_spatial_search.SearchInRadius(point_model_part.NodesBegin(), point_model_part.NodesEnd(), 0.3, results);

    // Print out results of the search
    for (auto& result : results) {
        std::cout << "Search result - Found: " << result.IsObjectFound() << ", Number of Global Results: " << result.NumberOfGlobalResults() << std::endl;
    }
}
```

#### Key points:

1. **Model setup**: This example creates a cubic skin model part using utility functions. In a real MPI run the `DataCommunicator` should come from the model part's `Communicator` (`rModelPart.GetCommunicator().GetDataCommunicator()`), and the model part must have been filled in parallel first (e.g. with `ParallelFillCommunicator`).

2. **Search wrapper initialization**: `ParallelSpatialSearch<GeometricalObjectsBins>` is initialized with the cube skin's elements and the data communicator. The elements to be searched are fixed at construction time.

3. **Search operation**: A search is performed for the point in a radius of 0.3 units. The search operation populates the results vector, which is then iterated to output the search findings.

### Python

For example if we have the following mesh divided in two partitions:

![](https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Search/search_wrapper_example.png?raw=true)

We can search the nodes in a radius of 0.5 from the nodes in X = 0.

![](https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Search/search_wrapper_example_2.png?raw=true)

We can perform that with the following script (see [`test_parallel_spatial_search.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/tests/test_parallel_spatial_search.py) for the full test):

~~~py
# Create search
search = KM.ParallelSpatialSearchKDTreeNode(model_part.Nodes, data_comm)

# Define radius
radius = 0.5

# Search only in the x=0 points
sub_model_part = model_part.CreateSubModelPart("SubModelPart")
for node in model_part.Nodes:
    if node.X == 0.0:
        sub_model_part.AddNode(node)

# Nodes array search
results = search.SearchInRadius(sub_model_part.Nodes, radius)

# Only in partitions were results are found
number_search_results = results.NumberOfSearchResults()
if number_search_results > 0:
    # GetResultIndices is a bulk operation of SpatialSearchResultContainerVector:
    # it returns, for every input point, the list of ids found
    all_ids = results.GetResultIndices()
    for i in range(number_search_results):
        node_results = results[i]
        global_id = node_results.GetGlobalIndex()
        ids = all_ids[i]
        number_of_global_results = node_results.NumberOfGlobalResults()
        if global_id > 0: # Solution defined in this rank
            rank = data_comm.Rank()
            print("Global id: ", global_id, " Rank: ", rank, " Number of local results: ", node_results.NumberOfLocalResults())
            if rank == 0:
               print("Global id: ", global_id, " Number of global results: ", number_of_global_results, " Ids: ", ids)
~~~

The output will look something like the following:

~~~sh
Global id:  1  Rank:  0  Number of local results:  2
Global id:  1  Rank:  1  Number of local results:  3
Global id:  1  Number of global results:  5  Ids:  [5, 4, 1, 2, 3]
Global id:  2  Rank:  0  Number of local results:  3
Global id:  2  Rank:  1  Number of local results:  4
Global id:  2  Number of global results:  7  Ids:  [5, 9, 4, 6, 1, 2, 3]
Global id:  6  Rank:  0  Number of local results:  4
Global id:  6  Rank:  1  Number of local results:  2
Global id:  6  Number of global results:  6  Ids:  [12, 5, 9, 4, 6, 2]
Global id:  12  Rank:  0  Number of local results:  5
Global id:  12  Rank:  1  Number of local results:  2
Global id:  12  Number of global results:  7  Ids:  [12, 5, 16, 9, 14, 6, 20]
Global id:  16  Rank:  0  Number of local results:  4
Global id:  16  Rank:  1  Number of local results:  1
Global id:  16  Number of global results:  5  Ids:  [12, 16, 9, 14, 20]
~~~

## Parameters & Defaults

The parameters and defaults, as returned by `GetDefaultParameters` and validated against the `Parameters` passed at construction, are the following:

~~~json
{
    "allocation_size"   : 1000,
    "bucket_size"       : 10
}
~~~

- `"allocation_size"`: An integer representing the initial allocation size used when collecting local `SearchInRadius` results per point. Default value is `1000`.
- `"bucket_size"`: An integer specifying the bucket size used when building the tree/bins-based search backends (not used by `GeometricalObjectsBins`). Default value is `10`.
