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
- **SYNCHRONOUS**: This mode presumably ensures a consistent, uniform handling of communication across all nodes or elements during spatial searches, regardless of their heterogeneity or distribution.
- **ASYNCHRONOUS**: Asynchronous communication, which is **not implemented yet**.

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

    // Initialize the results for these indexes
    container.InitializeResults(indexes.size());

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
bool is_local = true; // Flag to determine if the current operation is local to this processor (in parallel environments).
SpatialSearchResultContainer<GeometricalObject>* p_container = nullptr; // Pointer to a container for storing search results.

// The search results
SpatialSearchResultContainerVector<GeometricalObject> geometrical_object_results;

// Perform search, for example SearchNearest
... // Here the SearchNearest is performed, look at parallel_spatial_search

// Data communicator
ModelPart& r_root_model_part = rModelPart.GetRootModelPart();
const Communicator& r_communicator = r_root_model_part.GetCommunicator();
const DataCommunicator& r_data_communicator = r_communicator.GetDataCommunicator();

// Nodes array
auto& r_array_nodes = rModelPart.Nodes();
const auto it_node_begin = r_array_nodes.begin();

// Getting the number of results
const std::size_t number_of_search_results = geometrical_object_results.NumberOfSearchResults();

// Get all indices
std::vector<std::vector<std::vector<IndexType>>> all_indices;
geometrical_object_results.GetResultNodeIndices(all_indices);

// Get all is inside flags
std::vector<std::vector<bool>> all_is_inside;
geometrical_object_results.GetResultIsInside(all_is_inside, r_array_nodes, r_data_communicator, IsInsideTolerance);

// Get all shape functions
std::vector<std::vector<Vector>> all_shape_functions;
geometrical_object_results.GetResultShapeFunctions(all_shape_functions, r_array_nodes, r_data_communicator);

// Iterate overs solutions
for (unsigned int Index = 0; Index < number_of_search_results; ++Index) {
    // Get container
    p_container = &geometrical_object_results[Index];

    // If local point
    const bool is_local_point = p_container->IsLocalPoint();

    // At least one result
    if (p_container->NumberOfGlobalResults() > 0) {
        // Check
        KRATOS_DEBUG_ERROR_IF(p_container->NumberOfGlobalResults() > 1 && is_local_point) << "More than one result found for node " << (it_node_begin + p_container->GetLocalIndex())->Id() << " when using SearchNearest!" << std::endl;

        // If local search
        const bool is_local_search = p_container->IsLocalSearch(0);

        // Getting if it is local
        is_local = is_local_search && is_local_point;
        is_local = r_data_communicator.MaxAll(is_local); // At least active in one partition
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
                const bool is_found = all_is_inside[Index];

                // If found
                if (is_found) {
                    // Reserve
                    const std::size_t vector_size = rListVariables.size() * r_geometry.size();
                    partial_constraints.reserve(vector_size);
                    partial_constraints_pairs.reserve(vector_size);

                    // Compute shape function
                    const Vector& r_shape_function = all_shape_functions[Index][0];
                    for (unsigned int i = 0; i < r_geometry.size(); i++) {
                        if (std::abs(r_shape_function[i]) > ZeroTolerance) {
                            // Do something with the shape function
                        }
                    }
                }
            }
        } else {
            // Check if is inside
            const bool is_found = all_is_inside[Index];

            if (is_found) {
                // Getting the shape function and indices of the nodes, and parent index to be used linked to the results
                const Vector& r_shape_function = all_shape_functions[Index][0];
                const auto& r_indices = all_indices[Index][0];

                if (is_local_point) {
                    // Get the node
                    auto it_node = it_node_begin + p_container->GetLocalIndex();
                    auto p_node = *(it_node.base());

                    // Now that all nodes exist we can then generate the MPC
                    for (unsigned int i = 0; i < r_indices.size(); i++) {
                        const std::size_t index = r_indices[i];
                        if (std::abs(r_shape_function[i]) > ZeroTolerance) {
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
   - The appropriate search results container is accessed within the distributed computation loop.

##### Detecting Local Results

1. **Local vs global results**:
   - Each search result container (`SpatialSearchResultContainer`) contains methods like `IsLocalPoint()` and `NumberOfGlobalResults()` to determine if the current result is local to the MPI rank or if it's a result visible globally across multiple ranks.

2. **Handling local search**:
   - Local flags (`is_local_point`, `is_local_search`) are used to determine the execution flow and locality of data. This involves verifying node existence, geometry inclusion, and constraint application based on local data without needing further synchronization with other ranks.

3. **Using MPI**:
   - For results that might influence or be influenced by other ranks, MPI operations are used to synchronize and aggregate data (`r_data_communicator.MaxAll()` for coordinate comparison and aggregation).

4. **Distinction between local and remote data**:
   - Decisions based on local vs. remote data (e.g., whether a node should be considered isolated or constraints applied) hinge on the outcome of these local/global checks and MPI communications.

##### Practical steps for developers:

- **Initial setup**: Configure your environment, ensuring that the infrastructure for MPI communication is correctly set up.
- **Debugging and error checking**: Use assertions (like `KRATOS_DEBUG_ERROR_IF`) to ensure that results and operations make sense (e.g., no more than one result when expected, nodes are present in the model part).
- **Optimization**: Be mindful of performance implications, especially with MPI operations, data synchronization, and memory management (e.g., `shrink_to_fit`).

### Python

Similar to the previous C++ example, here you have a python example:

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestSpatialSearchResultContainerVector(KratosUnittest):
    def test_vector_workflow(self):
        data_comm = Kratos.Testing.GetDefaultDataCommunicator()
        rank = data_comm.Rank()
        world_size = data_comm.Size()
        num_searches_per_rank = 2

        model = Kratos.Model()
        model_part = model.CreateModelPart("Main")

        # 1. Create the container vector
        results_vector = Kratos.SpatialSearchResultContainerVectorGeometricalObject()

        # 2. Initialize containers for all local search points
        results_vector.InitializeResults(num_searches_per_rank)
        self.assertEqual(results_vector.NumberOfSearchResults(), num_searches_per_rank)

        # 3. Simulate a search and add results to each container
        # In a real case, these objects would be in a ModelPart
        for i in range(num_searches_per_rank):
            # Create a unique object ID based on rank and search index
            object_id = rank * num_searches_per_rank + i + 1
            node = model_part.CreateNewNode(object_id, 0, 0, 0)
            found_object = Kratos.GeometricalObject(Kratos.Point2D(node))

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
    * A `SpatialSearchResultContainerVector` is created. Using `InitializeResults(N)`, it's pre-sized to hold `N` empty `SpatialSearchResultContainer` instances, one for each search point this rank is responsible for.

2.  **Populating Local Results**
    * The script loops through the initialized containers. For each one, it simulates a search by creating a sample `GeometricalObject` and adding it to the corresponding container using `results_vector[i].AddResult(...)`. At this point, all data is still local to its rank.

3.  **Bulk Synchronization**
    * A single call to `SynchronizeAll(data_comm)` performs the main communication step. It efficiently gathers all results from all sub-containers across all MPI ranks and distributes the complete global picture back to everyone.

4.  **Bulk Data Retrieval**
    * After synchronization, a method like `GetDistances()` is called. This is the key advantage of the vector class: it fetches the distances for **all** search points and **all** their results in one optimized operation, returning a list of lists.