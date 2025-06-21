---
title: BinsDynamic
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: A dynamic binning data structure template for organizing and querying points in multi-dimensional space.
---

# Bins Dynamic

## Description

The `BinsDynamic` class template provides a dynamic binning data structure for organizing and querying points in a multi-dimensional space. It is parameterized by the dimension of the space, the point type, the container type for storing points, and other optional template parameters for specifying point iterators, distance iterators, and distance functions.

This class inherits from `TreeNode` to leverage common functionalities for tree-based data structures.

Also see [Bins (static)](bins_static) and [Geometrical Object Bins](geometrical_object_bins).

## Implementation

Can be found in [`bins_dynamic.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/bins_dynamic.h).

Derives from [`Tree`](tree.md) and [`TreeNode`](tree.md) class.

### Template arguments

- `TDimension`: The dimensionality of the space.
- `TPointType`: The type representing points in the space.
- `TContainerType`: The container type for storing points.
- `TPointerType`: The type of pointers or iterators to points in the container (default is inferred from `TContainerType`).
- `TIteratorType`: The type of iterators for traversing the container (default is inferred from `TContainerType`).
- `TDistanceIteratorType`: The type of iterators for storing distance values (default is `std::vector<double>::iterator`).
- `TDistanceFunction`: The type of the distance function used for querying (default is `Kratos::SearchUtils::SquaredDistanceFunction`).

### Search methods

- `SearchInRadius`: Iterates through potential cells based on the query point and radius, checking each object within these cells to gather those within the specified radius.
- `SearchInBox`: Iterates through potential cells based on the query point and bounding box points, checking each object within these cells to gather those within the specified bounding box.
- `SearchNearestPoint`: Focuses on finding the nearest object independently of the distance.

## Example usage

### C++

The following example of use for a simple case:

~~~c++
// Defining the point type for the search
using NodesContainerType = ModelPart::NodesContainerType;
using ResultNodesContainerType = NodesContainerType::ContainerType;
using VectorResultNodesContainerType = std::vector<ResultNodesContainerType>;
using PointType = PointObject<Node>;
using PointVector = std::vector<PointType::Pointer>;
using DistanceType = std::vector<double>;
using VectorDistanceType = std::vector<DistanceType>;

// Defining the PointVector
VectorResultNodesContainerType nodes_results;
VectorDistanceType result_distance;
PointVector points = SearchUtilities::PrepareSearch(rStructureNodes, rInputNodes, nodes_results, result_distance);

// Definitions
const int allocation_size = 1000;

/// BinsDynamic definitions
using DynamicBins = BinsDynamic<3ul, PointType, PointVector>;

// Creating the bins
DynamicBins dynamic_bins(points.begin(), points.end());

// Performing search
double radius = 1.0;
SearchUtilities::ParallelSearch(rInputNodes, radius, dynamic_bins, nodes_results, result_distance, allocation_size);
~~~
