---
title: Tree
keywords: search spatial_container core
tags: [search spatial_container]
sidebar: kratos_core_search
summary: A tree is a widely used abstract data type that represents a hierarchical tree structure with a set of connected nodes.
---

# Tree

![A generic, and so non-binary, unsorted, and with duplicated labels; arbitrary diagram of a tree](https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Search/tree.png?raw=true)

## Description

A tree is a widely used abstract data type that represents a hierarchical tree structure with a set of connected nodes. Each node in the tree can be connected to many children (depending on the type of tree), but must be connected to exactly one parent, except for the root node, which has no parent (i.e., the root node as the top-most node in the tree hierarchy).

## Implementation

Can be found in [`tree.h`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/spatial_containers/tree.h).

The implementation is mainly done in `Tree` and `TreeNode` classes.

### `Tree`

The `Tree` class uses the `TreeNode` as building blocks to create a more complex tree structure, such as binary trees, quad-trees, or oct-trees depending on the dimension and partition strategy specified by the template parameter `TPartitionType`.

#### Implementation details

- **Partitioning and Construction**: The constructor of the `Tree` class uses a partitioning strategy (`TPartitionType`) to organize points into a tree structure. It calculates a bounding box that encompasses all points and then recursively divides this space according to the strategy defined by `TPartitionType`.
- **Spatial Queries**: Implements methods to search for points within a specified radius, nearest points, and points within a bounding box, leveraging the capabilities of its `TreeNode` nodes.
- **Iterators**: Utilizes iterators to manage the collection of points, which allows the tree to handle large datasets efficiently.
- **Bounding Box**: Manages a bounding box around all points to optimize spatial queries, reducing the number of comparisons needed for a search.

#### Search methods

- `SearchInRadius`: Iterates through potential cells based on the query point and radius, checking each object within these cells to gather those within the specified radius.
- `SearchInBox`: Iterates through potential cells based on the query point and bounding box points, checking each object within these cells to gather those within the specified bounding box.
- `SearchNearestPoint`: Focuses on finding the nearest object independently of the distance.

### `TreeNode`

The `TreeNode` class is a generic node in a tree data structure that organizes spatial data through various geometric and partitioning methods. The class is highly customizable through its template parameters, which define the types of points, pointers, iterators, and distances used in the tree operations.

#### Template Parameters

- `TDimension`: Specifies the dimensionality of the space (e.g., 2D, 3D).
- `TPointType`: Represents the type of the points stored in the tree.
- `TPointerType`: Type of pointers to the objects stored.
- `TIteratorType`: Type of iterators used in the tree.
- `TDistanceIteratorType`: Type of iterators used for distances.
- `TIteratorIteratorType`: Typically an iterator for a vector of iterators.

#### Core Functionality

- Virtual functions like `SearchNearestPoint`, `SearchInRadius`, and `SearchInBox` are declared but not implemented in the `TreeNode`. These functions are intended to be defined in derived classes that specify particular tree behaviors like KD-tree, Octree, etc.
- This approach allows the `TreeNode` to be extremely flexible and adaptable to different kinds of spatial data structures by inheriting and implementing specific search algorithms.

## Usage

Examples in from [`KDTree`](kd_tree.md), [`OcTree`](octree.md), [`BinsDynamic`](bins_dynamic.md) and [`BinsStatic`](bins_static.md)  classes.
