---
title: GeometricalObjectsBins
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: A bins container for 3 dimensional GeometricalObject entities.
---

# Geometrical Object Bins

## Description

It provides efficient search in radius and search nearest methods. All of the geometries should be given at construction time. After constructing the bins the geometries cannot be modified. In case of any modification, the bins should be reconstructed.

Also see [Bins Dynamic](bins_dynamic) and [Bins (static)](bins_static).

## Implementation

Can be found in [`geometrical_object_bins.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/geometrical_object_bins.h).

### Constructors

Initialize the bins with a range of geometrical objects, automatically calculating the bounding box and appropriate cell sizes based on the objects' spatial distribution.

### Search methods

- `SearchInRadius`: Iterates through potential cells based on the query point and radius, checking each object within these cells to gather those within the specified radius.
- `SearchNearestInRadius`: Similar to `SearchInRadius` but focuses on finding the nearest object within the given radius, optimizing the search by incrementally increasing the search radius until the nearest object is found or all cells are covered.
- `SearchNearest`: Focuses on finding the nearest object independently of the distance.
- `SearchIsInside`: Determines if a point is directly inside any geometrical object, which can be useful for point-inclusion tests.

### Python exposition

1. **Constructors**: The class provides several overloaded constructors to initialize the `GeometricalObjectsBins` with either an `ElementsContainerType` or a `ConditionsContainerType`. Optionally also taking a `double` which represent a tolerance.

2. **Member functions**:

- `GetBoundingBox`: Exposed to retrieve the bounding box of the geometrical objects.
- `GetCellSizes`: Exposed to get the sizes of cells in the bins.
- `GetNumberOfCells`: Returns the number of cells.
- `GetTotalNumberOfCells`: Returns the total number of cells in all dimensions.

3. **Search functions**:

- `SearchInRadius`: Two overloaded versions are exposed. One takes a single point and a radius, performing a search within that radius and returning the results in a Python list. The other version accepts a container of nodes and a radius, returning a nested list where each sublist corresponds to the results for a particular node.
- `SearchNearestInRadius`: Similar to `SearchInRadius`, but presumably returns only the nearest object within the specified radius. It is also overloaded for single points and node containers.
- `SearchNearest`: Exposed for searching the nearest geometrical object to a point or each point in a container of nodes.
- `SearchIsInside`: Determines whether a point or nodes in a container are inside the geometric area defined by the geometrical objects.

Each search function captures lambda functions which interact with the `GeometricalObjectsBins` instance methods to perform searches, process results, and return them in a Python-friendly format using `py::list`.

## Example usage

### C++

The following is an example how to use it in C++ for example search all the results available in a given radius.

~~~c++
// Generate the cube skin
// Here, we create a skin model part of a cube using the utility function `CreateCubeSkinModelPart`.
// This function takes a reference to the current model and returns a reference to the newly created 'ModelPart'.
ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model);

// Create bins for spatial search
// This line initializes an object `bins` of type `GeometricalObjectsBins`. It is constructed using iterators 
// that define a range over the elements in 'r_skin_part'. This allows performing spatial searches on these elements.
GeometricalObjectsBins bins(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

// Define a vector to store search results
// A vector `results` is defined to store the results from the spatial search.
// The type stored in the vector is `GeometricalObjectsBins::ResultType`.
std::vector<GeometricalObjectsBins::ResultType> results;

// Define the center point for the search
// A `Point` object `center_point` is defined with coordinates (0.0, 0.0, 0.0).
Point center_point{0.0, 0.0, 0.0};

// Perform a radius search
// This function call performs a search within a radius of 0.29 units around `center_point`.
// The results of the search are stored in the `results` vector.
bins.SearchInRadius(center_point, 0.3, results);
~~~

### Python

Here an example how to use it in python:

```py
# Importing the Kratos Library
import KratosMultiphysics as KM

# Create a model
model = KM.Model()
# Create a model part
model_part = model.CreateModelPart("SampleModelPart")
model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3  # Define the spatial dimension

# Add variables that might be needed for the conditions
model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

# Create nodes
model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
model_part.CreateNewNode(3, 0.0, 1.0, 0.0)

# Create a triangular condition (assuming 2D domain for simplicity)
prop = model_part.GetProperties()[1]  # Default properties
model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1, 2, 3], prop)

# Example of using GeometricalObjectBins for searching conditions

# Initialize the search structure
search = KM.GeometricalObjectsBins(model_part.Conditions)

# Create a search node
search_node = KM.Node(4, 0.5, 0.1, 0.0)

# Search in radius
radius = 1.0
found_conditions = search.SearchInRadius(search_node, radius)

# Print found conditions
print("Conditions found in radius:", radius)
for cond in found_conditions:
    print("Condition ID:", cond.Id)
```
