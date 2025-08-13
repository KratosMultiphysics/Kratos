---
title: Sparse Linear Algebra (native Kratos implementation)
keywords: 
tags: [Sparse-Linear-Algebra-(native-Kratos-implementation).md]
sidebar: kratos_for_developers
summary: 
---

WORK IN PROGRESS

The Kratos prvides a basic implementation of sparse linear algebra. The implementation is designed to ease FEM operations and to provide a similar interface both in the SMP and MPI case.

# Construction and assembly of CSR Matrices (system matrices)
The library provides basic capabilities for the construction of CSR matrices by FEM assembly. 
A classical challenge for the construction of CSR matrices is that the graph of the matrix needs to be provided at once.
The kratos implementation provides a "graph" object designed to allow the sparse insertion of items. CSR matrices can be constructing by providing such a graph.

## SMP case

## MPI case

# Construction and assembly of Vectors (system vectors)
While the construction of Vectors ready for FEM assembly provides no difficulty in the SMP case, the effective use in a MPI context requires storing the communication patterns. The exact 


 
