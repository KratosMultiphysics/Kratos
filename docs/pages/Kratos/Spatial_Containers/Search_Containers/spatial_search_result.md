---
title: SpatialSearchResult
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: Class result of the spatial searches.
---

# Spatial Search Result

## Description

This class is the result of the spatial searches. It provides:
  - A global pointer to the object found.
  - Distance to the object if `IsDistanceCalculated()` is true
  - IsObjectFound if for example search nearest fails or not

This class is essential for applications involving spatial queries where it is crucial to not only find an object but also to know the distance to the object and handle cases where the search may fail. The use of `GlobalPointer` allows the system to work efficiently in a distributed environment, making it highly suitable for large-scale simulations handled by *Kratos*.

By abstracting the result of spatial searches into a dedicated class, *Kratos* provides a clean and reusable interface that enhances code maintainability and readability. This approach is beneficial in multi-physics simulations where different types of data and interactions need to be managed coherently.

See [Spatial Search Result Container](spatial_search_result_container) and [Spatial Search Result Container Vector](spatial_search_result_container_vector).

## Implementation

Can be found in [`spatial_search_result.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/spatial_search_result.h).

### Class: `SpatialSearchResult`
A template class designed to handle the results of spatial searches. The class is templated so it can handle to hold any kind of object.

#### Constructors
- **Default Constructor**: Initializes member variables to represent an unsuccessful search.
- **Parameterized Constructor**: Accepts a pointer to `TObjectType` and initializes based on whether the object exists. Default rank is zero, assuming a serial set up.

#### Operations
- **`Reset`**: Resets the search result to its initial state, indicating no object found and no distance calculated.

#### Accessors
- **Getters and Setters**: For accessing and modifying the object, its distance, and status flags (object found, distance calculated).

#### Inquiry
- **Status Checks**: Methods to check if the object was found and if the distance was calculated.

### Python exposition

Here’s a detailed explanation of how the C++ class `SpatialSearchResult` is exposed to Python:

#### Class specialization

Here's a breakdown of the objects that are specialized and their names as exposed in Python:

1. **Node** - This is specialized for spatial search and is exposed in Python as `"SpatialSearchResultNode"`.
2. **GeometricalObject** - This is another specialization for spatial search results and is exposed in Python as `"SpatialSearchResultGeometricalObject"`.
3. **Element** - Specialized for spatial search results concerning elements, and it is exposed in Python as `"SpatialSearchResultElement"`.
4. **Condition** - This type is also specialized for spatial search results, specifically for conditions, and is exposed in Python as `"SpatialSearchResultCondition"`.

#### Member functions
 - `Reset`: Binds the `Reset` method that resets the internal state of the search result.
 - `Get`: Binds the `Get` method, which returns the object stored in the result.
 - `Set`: Binds the `Set` method to set the object in the result.
 - `GetDistance`: Binds the `GetDistance` method to retrieve the distance at which the object was found.
 - `SetDistance`: Binds the `SetDistance` method to set the distance.
 - `IsObjectFound`: Binds a method to check if an object was found in the search.
 - `IsDistanceCalculated`: Binds a method to check if the distance has been calculated.

## Example usage

### Python

~~~python
import KratosMultiphysics as KM

# Create a new SpatialSearchResult for a specific type, say Node
result = KM.SpatialSearchResultNode()

# Use the methods bound to the SpatialSearchResultNode class
result.Set(node_pointer)  # Set a node
print(result.Get())  # Get the node
print(result.GetDistance())  # Get the distance
result.SetDistance(10.0)  # Set the distance
print(result.IsObjectFound())  # Check if the object was found
print(result.IsDistanceCalculated())  # Check if the distance was calculated
~~~