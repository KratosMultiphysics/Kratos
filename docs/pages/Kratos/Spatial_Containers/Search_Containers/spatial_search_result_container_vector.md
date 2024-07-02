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

#### Instantiation

The class is exposed with several template instantiations for different types (`Node`, `GeometricalObject`, `Element`, `Condition`) with two communication types:
- **SYNCHRONOUS_HOMOGENEOUS**: This mode presumably ensures a consistent, uniform handling of communication across all nodes or elements during spatial searches, regardless of their heterogeneity or distribution.
- **SYNCHRONOUS_HETEROGENEOUS**: In contrast, this mode might allow for different handling or processing strategies depending on the node or element type or state, providing more flexibility in operations that require acknowledgment of diversity in the dataset or environment.

Each type and communication mode combination is then bound to the module with a specific class name, such as `SpatialSearchResultContainerVectorNode` for nodes with homogeneous communication, or `SpatialSearchResultContainerVectorNodeHeterogeneous` for nodes with heterogeneous communication.

#### Methods exposed

##### Constructor (`__init__`)
Initializes a new instance of `SpatialSearchResultContainerVector`.

##### NumberOfSearchResults
Returns the number of search results contained in the vector.

##### InitializeResult
Initializes the search results to a default state.

##### HasResult
Checks if a result is available in the container.

##### Clear
Clears all the contents of the container, resetting it to its initial state.

##### SynchronizeAll
Synchronizes all search results across different computational entities or nodes, ensuring consistency. In order to be efficient, and this is a very important operation that would affect performance, all the global pointers of all the solutions stored are processed at once as it would be stored in this class.

See [Spatial Search Result Container SynchronizeAll](spatial_search_result_container#SynchronizeAll).

##### GetDistances
Retrieves the distances associated with the search results. This method can be called instead of the one in [Spatial Search Result Container](spatial_search_result_container) in order to perform all the global pointers operations at once for the sake of efficiency.

See [Spatial Search Result Container GetDistances](spatial_search_result_container#GetDistances).

##### GetResultIsLocal
Checks if the result is local to the computational entity or node. This method can be called instead of the one in [Spatial Search Result Container](spatial_search_result_container) in order to perform all the global pointers operations at once for the sake of efficiency.

See [Spatial Search Result Container GetResultIsLocal](spatial_search_result_container#GetResultIsLocal).

##### GetResultRank
Retrieves the rank of the result within the container. This method can be called instead of the one in [Spatial Search Result Container](spatial_search_result_container) in order to perform all the global pointers operations at once for the sake of efficiency.

See [Spatial Search Result Container GetResultRank](spatial_search_result_container#GetResultRank).

##### GetResultIsActive
Determines whether the result is active. This method can be called instead of the one in [Spatial Search Result Container](spatial_search_result_container) in order to perform all the global pointers operations at once for the sake of efficiency.

See [Spatial Search Result Container GetResultIsActive](spatial_search_result_container#GetResultIsActive).

##### GetResultIsInside
Checks if the result is inside a specified boundary or criteria. This method can be called instead of the one in [Spatial Search Result Container](spatial_search_result_container) in order to perform all the global pointers operations at once for the sake of efficiency.

See [Spatial Search Result Container GetResultIsInside](spatial_search_result_container#GetResultIsInside).

##### GetResultShapeFunctions
Retrieves the shape functions associated with the search result. This method can be called instead of the one in [Spatial Search Result Container](spatial_search_result_container) in order to perform all the global pointers operations at once for the sake of efficiency.

See [Spatial Search Result Container GetResultShapeFunctions](spatial_search_result_container#GetResultShapeFunctions).

##### GetResultIndices
Gets the indices of the search results. This method can be called instead of the one in [Spatial Search Result Container](spatial_search_result_container) in order to perform all the global pointers operations at once for the sake of efficiency.

See [Spatial Search Result Container GetResultIndices](spatial_search_result_container#GetResultIndices).

##### GetResultNodeIndices
Retrieves the node indices of the search results. This method can be called instead of the one in [Spatial Search Result Container](spatial_search_result_container) in order to perform all the global pointers operations at once for the sake of efficiency.

See [Spatial Search Result Container GetResultNodeIndices](spatial_search_result_container#GetResultNodeIndices).

##### GetResultPartitionIndices
Gets the partition indices related to the search results. This method can be called instead of the one in [Spatial Search Result Container](spatial_search_result_container) in order to perform all the global pointers operations at once for the sake of efficiency.

See [Spatial Search Result Container GetResultPartitionIndices](spatial_search_result_container#GetResultPartitionIndices).

##### GetResultCoordinates
Retrieves the coordinates of the search results. This method can be called instead of the one in [Spatial Search Result Container](spatial_search_result_container) in order to perform all the global pointers operations at once for the sake of efficiency.

See [Spatial Search Result Container GetResultCoordinates](spatial_search_result_container#GetResultCoordinates).

##### GetGlobalPointerCommunicator
Obtains the communicator used for global pointer operations across computational entities.

See [Spatial Search Result Container GetGlobalPointerCommunicator](spatial_search_result_container#GetGlobalPointerCommunicator).

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
#include "geometries/point.h"
#include "spatial_containers/spatial_search_result_container_vector.h"

int main() {
    // Assuming you have a suitable MPI environment set up
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    // Create an instance of the container for GeometricalObjects
    SpatialSearchResultContainerVector<GeometricalObject> container;

    // Assume 'indexes' are the indices for which you want to initialize results
    std::vector<std::size_t> indexes = {0, 1, 2};  // Example indexes
    std::vector<const DataCommunicator*> data_comms(indexes.size(), &r_data_comm);

    // Initialize the results for these indexes
    container.InitializeResults(data_comms);

    // Now you can add results to these initialized positions
    for (std::size_t idx : indexes) {
        GeometricalObject* pObj = new GeometricalObject;  // Assuming GeometricalObject is properly defined
        SpatialSearchResult<GeometricalObject> result(pObj);
        result.SetDistance(0.1 * idx);  // Example distance setting

        container[idx].AddResult(result);  // Adding result to the specific container
    }

    // Optionally, synchronize results across MPI processes
    container.SynchronizeAll(r_data_comm);

    // Access and print out results
    for (std::size_t idx : indexes) {
        auto& results = container[idx].GetLocalResults();
        for (auto& res : results) {
            std::cout << "Object at index " << idx << " has distance: " << res.GetDistance() << std::endl;
        }
    }

    // Remember to clean up any dynamic allocations if needed
    for (std::size_t idx : indexes) {
        auto& results = container[idx].GetLocalResults();
        for (auto& res : results) {
            delete res.GetResult();  // Assuming GetResult returns a pointer to GeometricalObject
        }
    }

    return 0;
}
```

### Key Points to Remember
- You need an MPI environment to use the `SpatialSearchResultContainerVector` as it handles results across different processes.
- You should initialize the results for the specific indexes you are going to use.
- `AddResult` method allows you to add spatial search results into the container.
- You can synchronize results across different processes using the `SynchronizeAll` method if needed.
- Always manage memory carefully, especially when dynamically allocating objects like `GeometricalObject`.

#### Advanced example

Here a code example that primarily focuses on performing spatial searches to identify and process geometrical objects relative to specified nodes within a model part. It integrates with search algorithms, data communication across computational partitions, and error-checking mechanisms.

```cpp
// Auxiliary stuff
Vector shape_function; // Stores shape functions used in finite element analysis.
Point::CoordinatesArrayType aux_coords; // Auxiliary coordinates for internal computations.
array_1d<double, 3> point; // Represents a 3D point in space.
bool is_local = true; // Flag to determine if the current operation is local to this processor (in parallel environments).
SpatialSearchResultContainer<GeometricalObject>* p_container = nullptr; // Pointer to a container for storing search results.

// The search results
SpatialSearchResultContainerVector<GeometricalObject> geometrical_object_results;

// Perform search, for example SearchNearest
... // Here the SearchNearest is performed, look at search_wrapper

// Nodes array
auto& r_array_nodes = rModelPart.Nodes();
const auto it_node_begin = r_array_nodes.begin();

// Getting the number of results
const std::size_t number_of_search_results = geometrical_object_results.NumberOfSearchResults();

// Iterate overs solutions
for (unsigned int Index = 0; Index < number_of_search_results; ++Index) {
    // Get container
    p_container = &geometrical_object_results[Index];

    // Using barrier method to avoid communication issues
    p_container->Barrier();

    // If local point
    const bool is_local_point = p_container->IsLocalPoint();

    // At least one result
    if (p_container->NumberOfGlobalResults() > 0) {
        // Check
        KRATOS_DEBUG_ERROR_IF(p_container->NumberOfGlobalResults() > 1 && is_local_point) << "More than one result found for node " << (it_node_begin + p_container->GetLocalIndex())->Id() << " when using SearchNearest!" << std::endl;

        // If local search
        const bool is_local_search = p_container->IsLocalSearch(0);

        // Get data communicator
        const auto& r_sub_data_communicator = p_container->GetDataCommunicator();

        // Getting if it is local
        is_local = is_local_search && is_local_point;
        is_local = r_sub_data_communicator.MaxAll(is_local); // At least active in one partition
        if (is_local) {
            if (is_local_point) {
                // Get the node
                auto it_node = it_node_begin + p_container->GetLocalIndex();
                auto p_node = *(it_node.base());
                auto& r_node = *p_node;

                // Check is actually inside (not just found)
                KRATOS_ERROR_IF(p_container->NumberOfLocalResults() == 0) << "Local results should not be null" << std::endl;
                const auto& r_local_result = p_container->GetLocalResults()[0];
                const GeometricalObject* p_geometrical_object = r_local_result.Get().get(); // Intrussive_ptr?
                const auto& r_geometry = p_geometrical_object->GetGeometry();
                const bool is_found = r_geometry.IsInside(r_node.Coordinates(), aux_coords, IsInsideTolerance);

                // If found
                if (is_found) {
                    // Reserve
                    const std::size_t vector_size = rListVariables.size() * r_geometry.size();
                    partial_constraints.reserve(vector_size);
                    partial_constraints_pairs.reserve(vector_size);

                    // Compute shape function
                    shape_function = r_geometry.ShapeFunctionsValues(shape_function, aux_coords);
                    for (unsigned int i = 0; i < r_geometry.size(); i++) {
                        if (std::abs(shape_function[i]) > ZeroTolerance) {
                            // Do something with the shape function
                        }
                    }
                }
            }
        } else {
            if (is_local_point) {
                auto it_node = it_node_begin + p_container->GetLocalIndex();
                noalias(point) = it_node->Coordinates();
            } else {
                for (unsigned int i = 0; i < 3; i++) {
                    point[i] = -MaxValue;
                }
            }
            for (unsigned int i = 0; i < 3; i++) {
                point[i] = r_sub_data_communicator.MaxAll(point[i]);
            }

            // Check if is inside
            const bool is_found = p_container->GetResultIsInside(point, IsInsideTolerance)[0];

            if (is_found) {
                // Getting the shape function and indices of the nodes, and parent index to be used linked to the results
                const auto shape_function = p_container->GetResultShapeFunctions(point)[0];
                const auto indices = p_container->GetResultNodeIndices()[0];

                if (is_local_point) {
                    // Get the node
                    auto it_node = it_node_begin + p_container->GetLocalIndex();
                    auto p_node = *(it_node.base());

                    // Now that all nodes exist we can then generate the MPC
                    for (unsigned int i = 0; i < indices.size(); i++) {
                        const std::size_t index = indices[i];
                        if (std::abs(shape_function[i]) > ZeroTolerance) {
                            KRATOS_DEBUG_ERROR_IF_NOT(rPrimaryModelPart.HasNode(index)) << "Node " << index << " has not being bring to the modelpart" << std::endl;
                            auto& r_node_linked = rPrimaryModelPart.GetNode(index);
                            // Do something with the shape function
                        }
                    }
                }
            }
        }
    }
}
```

##### How to use `geometrical_object_results`

1. **Retrieve search results**:
   - The number of search results is obtained using `geometrical_object_results.NumberOfSearchResults()`. This number defines the loop bounds for processing each search result.

2. **Access individual results**:
   - Depending on the type of search communication (`SYNCHRONOUS_HOMOGENEOUS` or `SYNCHRONOUS_HETEROGENEOUS`), the appropriate search results container is accessed within the distributed computation loop.

3. **Barrier synchronization**:
   - A barrier method is called (`rTLS.p_container->Barrier()`) to synchronize across processes, ensuring all necessary data is ready before proceeding, thus preventing premature access to partial data. This method is only actually used for `SYNCHRONOUS_HETEROGENEOUS`, but if you want your code to work in both `SYNCHRONOUS_HOMOGENEOUS` and `SYNCHRONOUS_HETEROGENEOUS` then should be called in all cases.

##### Detecting Local Results

1. **Local vs global results**:
   - Each search result container (`SpatialSearchResultContainer`) contains methods like `IsLocalPoint()` and `NumberOfGlobalResults()` to determine if the current result is local to the MPI rank or if it's a result visible globally across multiple ranks.

2. **Handling local search**:
   - Local flags (`is_local_point`, `is_local_search`) are used to determine the execution flow and locality of data. This involves verifying node existence, geometry inclusion, and constraint application based on local data without needing further synchronization with other ranks.

3. **Using MPI**:
   - For results that might influence or be influenced by other ranks, MPI operations are used to synchronize and aggregate data (`r_sub_data_communicator.MaxAll()` for coordinate comparison and aggregation).

4. **Distinction between local and remote data**:
   - Decisions based on local vs. remote data (e.g., whether a node should be considered isolated or constraints applied) hinge on the outcome of these local/global checks and MPI communications.

##### Practical steps for developers:

- **Initial setup**: Configure your environment to handle both types of searches (`SYNCHRONOUS_HOMOGENEOUS` and `SYNCHRONOUS_HETEROGENEOUS`), ensuring that the infrastructure for MPI communication is correctly set up.
- **Debugging and error checking**: Use assertions (like `KRATOS_DEBUG_ERROR_IF`) to ensure that results and operations make sense (e.g., no more than one result when expected, nodes are present in the model part).
- **Optimization**: Be mindful of performance implications, especially with MPI operations, data synchronization, and memory management (e.g., `shrink_to_fit`).
