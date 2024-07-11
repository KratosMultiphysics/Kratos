---
title: Bins
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: Class for searching and organizing spatial data using a binning technique.
---

# Bins (static)

## Description

The bins is deigned for searching and organizing spatial data using a binning technique. This technique refers to an approach where data is first grouped into bins or categories before the search is conducted. This method can be particularly useful for improving the efficiency of search processes in large datasets. The method is based on:

1. **Data segregation**: The data set is divided into multiple bins based on certain criteria, often the range of data values. Each bin contains data that falls within a specific interval.

2. **Efficient searching**: When a search query is issued, instead of searching the entire dataset, the algorithm first identifies the appropriate bin or bins where the search keys are likely to be found.

Also see [Bins Dynamic](bins_dynamic) and [Geometrical Object Bins](geometrical_object_bins).

## Implementation

Can be found in [`bins_static.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/bins_static.h).

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
const int bucket_size = 1000;

/// StaticBins definitions
using StaticBinsTree = Tree<Bins<3ul, PointType, PointVector>>;

// Creating the tree
StaticBinsTree static_bins_tree(points.begin(), points.end(), bucket_size);

// Performing search
double radius = 1.0;
SearchUtilities::ParallelSearch(rInputNodes, radius, static_bins_tree, nodes_results, result_distance, allocation_size);
~~~
