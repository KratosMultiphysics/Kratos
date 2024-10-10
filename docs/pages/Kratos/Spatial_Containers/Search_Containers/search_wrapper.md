---
title: SearchWrapper
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: This is a search wrapper ready for MPI searches
---

# Search Wrapper

## Description

This is a search wrapper, which must be adapted and specialized for every search object.The class is designed to efficiently handle searches across different processes in a distributed computing environment using MPI.

Spatial container are used to manage and query spatial data efficiently. Classes like [`GeometricalObjectsBins`](geometrical_object_bins) and various tree-based structures (e.g., [KDTree](kd_tree), [OCTree](octree), [BinsDynamic](bind_dynamic)) are used depending on the application requirements.

See [`SpatialSearchResultContainer`](spatial_search_result_container) and [`SpatialSearchResultContainerVector`](spatial_search_result_container_vector).

## Implementation

Can be found in [`search_wrapper.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/search_wrapper.h).

### Template arguments

- `TSearchObject`: The type of the spatial search object.
- `TSpatialSearchCommunication`: The communication strategy, which can be synchronous or asynchronous, and homogeneous or heterogeneous. Considers the enum `SpatialSearchCommunication`, see [`SpatialSearchResultContainer`](spatial_search_result_container):
  - `SYNCHRONOUS_HOMOGENEOUS`: All partitions are fully aware of the entire dataset.
  - `SYNCHRONOUS_HETEROGENEOUS`: Partitions are aware of sub-sets of the dataset.
  - `ASYNCHRONOUS`: Asynchronous communication, which is not implemented yet.

### Result containers

- **[`SpatialSearchResultContainer`](spatial_search_result_container) and [`SpatialSearchResultContainerVector`](spatial_search_result_container_vector)**: These containers store the results of spatial searches, handling different data types and ensuring compatibility with MPI.

### Python exposition

#### Template parameters:

- `TSearchObject`: Type of the search object.
- `TSpatialSearchCommunication`: Type of spatial search communication.

#### Exposed methods:

1. **Constructors**:
   - Different constructors are exposed based on the type of `ObjectType` within `TSearchObject`. Depending on whether the object is a `Node`, `Element`, or `Condition`, appropriate container types are used as arguments along with a `DataCommunicator` and optionally, `Parameters`.

2. **Spatial Search Methods**:
   - `GetGlobalBoundingBox`: Retrieves the global bounding box of the search space.
   - `SearchInRadius`: Performs a radius-based search.
   - `SearchNearestInRadius`: Searches for the nearest objects within a specified radius.
   - `SearchNearest`: Searches for the nearest objects.
   - `SearchIsInside`: Checks if objects are inside a specified space (only if `TSearchObject` is `GeometricalObjectBins`).

#### Instantiation

##### 1. `SearchWrapperGeometricalObjectBins`
   - For both `SYNCHRONOUS_HOMOGENEOUS` and `SYNCHRONOUS_HETEROGENEOUS` communication types.
   - Exposes constructors for `GeometricalObjectsBins` and the spatial search methods as described above.

##### 2. KDTree Wrappers
   - Different variants of KDTree for `Node`, `Element`, and `Condition`:
     - `SearchWrapperKDTreeNode`
     - `SearchWrapperKDTreeElement`
     - `SearchWrapperKDTreeCondition`
   - Both `SYNCHRONOUS_HOMOGENEOUS` and `SYNCHRONOUS_HETEROGENEOUS` versions.

##### 3. OCTree Wrappers
   - Analogous to KDTree wrappers, tailored for OCTree partitions:
     - `SearchWrapperOCTreeNode`
     - `SearchWrapperOCTreeElement`
     - `SearchWrapperOCTreeCondition`
   - Both `SYNCHRONOUS_HOMOGENEOUS` and `SYNCHRONOUS_HETEROGENEOUS` versions.

##### 4. Static and Dynamic Bins Wrappers
   - Exposes wrappers for static and dynamic bins structures for `Node`, `Element`, and `Condition`:
     - `SearchWrapperStaticBinsTreeNode`
     - `SearchWrapperDynamicBinsNode`, etc.
   - Both `SYNCHRONOUS_HOMOGENEOUS` and `SYNCHRONOUS_HETEROGENEOUS` versions.

## Example usage

### C++

Here's an example of how you can utilize the `SearchWrapperGeometricalObjectsBins` class in C++:

```cpp
#include "containers/model.h"
#include "tests/test_utilities/cpp_tests_utilities.h"
#include "spatial_containers/search_wrapper.h"
#include "mpi/utilities/parallel_fill_communicator.h"

using namespace Kratos;

static void Example() {
    Model current_model;

    // Create a cube model part
    ModelPart& cube_skin = CppTestsUtilities::CreateCubeSkinModelPart(current_model, 0.6, 0.9, 0.3);

    // Setup the parallel environment for the model part
    const DataCommunicator& data_communicator = Testing::GetDefaultDataCommunicator();
    ParallelFillCommunicator(cube_skin, data_communicator).Execute();

    // Define the type of the search wrapper using homogeneous communication
    using SearchWrapperType = SearchWrapperGeometricalObjectsBins;

    // Instantiate the search wrapper with the elements of the model part
    SearchWrapperType search_wrapper(cube_skin.Elements(), data_communicator);

    // Create a model part for points to be searched
    ModelPart& point_model_part = current_model.CreateModelPart("PointModelPart");
    point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // Adding a point in the point model part
    auto p_node = point_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;

    // Create a vector to store the search results
    typename SearchWrapperType::ResultContainerVectorType results;

    // Conduct a search in radius of 0.3 units around the newly added node
    search_wrapper.SearchInRadius(point_model_part.NodesBegin(), point_model_part.NodesEnd(), 0.3, results);

    // Print out results of the search
    for (auto& result : results) {
        std::cout << "Search result - Found: " << result.IsObjectFound() << ", Number of Global Results: " << result.NumberOfGlobalResults() << std::endl;
    }
}
```

#### Key points:

1. **Model setup**: This example creates a cubic skin model part using utility functions. The model dimensions are specified, and communication details are set up using MPI utilities, which is typical in distributed environments.

2. **Search wrapper initialization**: The `SearchWrapperGeometricalObjectsBins` is initialized with the cube skin's elements and the default data communicator.

3. **Search operation**: A search is performed for points in a radius of 0.3 units. The search operation populates the results vector, which is then iterated to output the search findings.

4. **MPI environment**: Depending on your setup, you might need to initialize and finalize the MPI environment to handle parallel computations properly. This is suggested in the comments but commented out for simplicity and context dependency.

### Python

For example if we have the following mesh divided in two partitions:

![](https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Search/search_wrapper_example.png?raw=true)

We can search the nodes in a radius of 0.5 from the nodes in X = 0.

![](https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Search/search_wrapper_example_2.png?raw=true)

We can perform that with the following script:

~~~py
# Create search
search = KM.SearchWrapperKDTreeNode(model_part.Nodes, data_comm)

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
    for i in range(5):
        node_results = results[i]
        global_id = node_results.GetGlobalIndex()
        ids = node_results.GetResultIndices()
        number_of_global_results = node_results.NumberOfGlobalResults()
        if global_id > 0: # Solution defined in this rank
            rank = data_comm.Rank()
            print("Global id: ", global_id, " Rank: ", rank, " Number of local results: ", node_results.NumberOfLocalResults())
            if rank == 0:
               print("Global id: ", global_id, " Number of global results: ", number_of_global_results, " Ids: ", ids)
~~~

The output will looks something like the following:

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

The parameters and defaults are the followings:

~~~json
{
    "allocation_size"   : 1000,
    "bucket_size"       : 10
}
~~~

- `"allocation_size"`: An integer representing the initial allocation size. Default value is `1000`.
- `"bucket_size"`: An integer specifying the size of each bucket. Default value is `10`.