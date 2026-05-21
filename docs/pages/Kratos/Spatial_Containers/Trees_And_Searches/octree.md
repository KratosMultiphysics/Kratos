---
title: OcTree
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: A tree where each node has up to eight children.
---

# OcTree

![Schematic drawing of an octree](https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Search/octree.png?raw=true)

## Description

An octree is essentially a tree where each node has up to eight children, hence the prefix "oct-". This structure is particularly useful for subdividing a 3D space into smaller segments. The root of an octree represents a cubic volume, and each of its eight children represents one of the eight subdivisions of that cube. This subdivision continues recursively, depending on the depth of the tree or until a certain condition is met (like a minimal volume size or a specific number of contained elements).

### Advantages of octrees

- **Efficient Storage**: By dividing space into hierarchical segments, octrees can store 3D objects more efficiently, reducing the complexity of intersection tests and other computations.
- **Scalability**: Octrees can handle very large amounts of data by adjusting the depth of the tree.
- **Fast Query Performance**: They offer efficient query performance for range searches, nearest neighbor searches, and collision detection, among other applications.

## Implementation

Can be found in [`octree.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/octree.h).

The implementation of an octree typically involves recursive subdivision of space, and storing objects in the nodes that completely contain them. Optimizations might include limiting the depth of the tree to prevent excessive recursion and balancing the tree to ensure that data is evenly distributed across nodes.

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

/// Octree definitions
using Octree = Tree<OCTreePartition<Bucket<3ul, PointType, PointVector>>>;

// Creating the tree
Octree octree(points.begin(), points.end(), bucket_size);

// Performing search
double radius = 1.0;
SearchUtilities::ParallelSearch(rInputNodes, radius, octree, nodes_results, result_distance, allocation_size);
~~~
