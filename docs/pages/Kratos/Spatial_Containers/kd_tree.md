---
title: KDTree
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: A k-d tree (short for k-dimensional tree) is a space-partitioning data structure for organizing points in a k-dimensional space.
---

# KDTree

![A 3-d tree](images/3dtree.png)

## Description

A k-d tree (short for k-dimensional tree) is a space-partitioning data structure for organizing points in a k-dimensional space. k-d trees are useful in searching for points in a spatial context, and nearest neighbor searches.

1. **Recursive binary tree structure**: A k-d tree is a binary tree in which every node represents an axis-aligned hyperrectangle in space. Each non-leaf node generates a splitting hyperplane that divides the space into two parts, known as half-spaces. Points to the left of this hyperplane are represented by the left subtree of that node and points to the right by the right subtree.

2. **Splitting dimensions and values**: The tree is constructed by recursively choosing one of the k dimensions of the space to use as the splitting dimension. Typically, the dimension is chosen in a round-robin fashion (cycling through dimensions) or based on the dimension with the greatest variance among the points yet to be stored. At each step, a point in the current set is selected to act as the pivot (often using the median value of the chosen dimension) to partition the dataset into two subsets.

3. **Construction**: Starting with the entire set of points, a dimension is chosen for splitting and a median is selected based on that dimension. This point becomes a node in the k-d tree, and the process is repeated recursively for the subset of points on each side of the median point until the subsets are empty, resulting in the creation of leaf nodes.

4. **Search operations**: Searching for a point in a k-d tree involves traversing the tree from the root, making decisions at each node based on how the point to be found compares to the node in the chosen dimension. This is efficient because it eliminates half of the search space with each level of tree depth.

5. **Nearest neighbour search**: This is one of the most common uses of k-d trees. The algorithm involves searching down the tree to find the leaf node closest to the target point, then unwinding the recursion, checking at each node whether there could be a closer point in the opposite half-space (which might intersect the hyperrectangle containing the target point).

6. **Complexity**: The average time complexity for constructing a k-d tree is \(O(n \log n)\), where \(n\) is the number of points, assuming the median is found efficiently. The search and insertion operations generally take \(O(\log n)\) time in well-balanced trees.

## Implementation

Can be found in [`kd_tree.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/kd_tree.h).

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

/// KDtree definitions
using KDTree = Tree<KDTreePartition<Bucket<3ul, PointType, PointVector>>>;

// Creating the tree
KDTree kd_tree(points.begin(), points.end(), bucket_size);

// Performing search
double radius = 1.0;
SearchUtilities::ParallelSearch(rInputNodes, radius, kd_tree, nodes_results, result_distance, allocation_size);
~~~
