---
title: Overview
keywords: 
tags: [Overview.md]
sidebar: statistics_application
summary: 
---
Statistics application consist of widely used methods to calculate statistics in various containers of KratosMultiphysics. There are mainly two groups of statistical methods namely Spatial and Temporal. Spatial methods calculate statistics on spatial containers and output the values whenever they are called. Temporal methods calculate statistics on the fly in a transient simulation. All the temporal methods gurantee that, the resultant statistics will be same as the ones if one calculates accumulating all the data upto that time instance and calculate the same statistics. All of these methods in each group is OpenMP and MPI compatible, and tested.

Following table summarize capabilities of Statistics Application.

| Statistical Methods                   | Norm Types              | Spatial Domain                                                | Temporal Domain                                                              | Data types |
|---------------------------------------|-------------------------|---------------------------------------------------------------|------------------------------------------------------------------------------|------------|
| [Sum](#sum)                           | [Value](#value)         | [Spatial methods](#spatial-methods)                           | [Temporal methods](#temporal-methods)                                        | Double     |
| [Mean](#mean)                         | [Magnitude](#magnitude) | [Spatial containers](#spatial-containers)                     | [Temporal containers](#temporal-containers)                                  | Array 3D   |
| [Root mean square](#root-mean-square) | [Euclidean](#euclidean) | [nodal_historical](#spatial-nodal-historical)                 | [nodal_historical_historical](#temporal-nodal-historical-historical)         | Vector     |
| [Variance](#variance)                 | [Infinity](#infinity)   | [nodal_non_historical](#spatial-nodal-non-historical)         | [nodal_historical_non_historical](#temporal-nodal-historical-non-historical) | Matrix     |
| [Min](#min)                           | [P-Norm](#p-norm)       | [condition_non_historical](#spatial-condition-non-historical) | [nodal_non_historical](#temporal-nodal-non-historical)                       |            |
| [Max](#max)                           | [Lpq-Norm](#lpq-norm)   | [element_non_historical](#spatial-element-non-historical)     | [element_non_historical](#temporal-element-non-historical)                   |            |
| [Median](#median)                     | [Frobenius](#frobenius) | [Spatial statistics process](#spatial-statistics-process)     | [condition_non_historical](#temporal-condition-non-historical)               |            |
| [Distribution](#distribution)         | [Trace](#trace)         |                                                               | [Temporal statistics process](#temporal-statistics-process)                  |            |
| [Norm methods](#norm-methods)         | [Index](#index-based)         |                                                               |                                                                              |            |
|                                       | [Component](#component-based) |                                                               |                                                                              |            |
