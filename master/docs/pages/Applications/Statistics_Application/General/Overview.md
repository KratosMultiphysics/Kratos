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
| [Sum](../Statistical_Methods/Sum.html)                           | [Value](../Norm_Methods/value.html)         | [Spatial methods](../Domains/Spatial/spatial_methods.html)                           | [Temporal methods](../Domains/Temporal/temporal_methods.html)                                        | Double     |
| [Mean](../Statistical_Methods/Mean.html)                         | [Magnitude](../Norm_Methods/magnitude.html) | [Spatial containers](../Domains/Spatial/spatial_containers.html)                     | [Temporal containers](../Domains/Temporal/temporal_containers.html)                                  | Array 3D   |
| [Root mean square](../Statistical_Methods/Root_Mean_Square.html) | [Euclidean](../Norm_Methods/euclidean.html) | [nodal_historical](../Domains/Spatial/spatial_containers.html#spatial-nodal-historical)                 | [nodal_historical_historical](../Domains/Temporal/temporal_containers.html#temporal-nodal-historical-historical)         | Vector     |
| [Variance](../Statistical_Methods/Variance.html)                 | [Infinity](../Norm_Methods/infinity.html)   | [nodal_non_historical](../Domains/Spatial/spatial_containers.html#spatial-nodal-non-historical)         | [nodal_historical_non_historical](../Domains/Temporal/temporal_containers.html#temporal-nodal-historical-non-historical) | Matrix     |
| [Min](../Statistical_Methods/Min.html)                           | [P-Norm](../Norm_Methods/p_norm.html)       | [condition_non_historical](../Domains/Spatial/spatial_containers.html#spatial-condition-non-historical) | [nodal_non_historical](../Domains/Temporal/temporal_containers.html#temporal-nodal-non-historical)                       |            |
| [Max](../Statistical_Methods/Max.html)                           | [Lpq-Norm](../Norm_Methods/lpq_norm.html)   | [element_non_historical](../Domains/Spatial/spatial_containers.html#spatial-element-non-historical)     | [element_non_historical](../Domains/Temporal/temporal_containers.html#temporal-element-non-historical)                   |            |
| [Median](../Statistical_Methods/Median.html)                     | [Frobenius](../Norm_Methods/frobenius.html) | [Spatial statistics process](../Domains/Spatial/spatial_statistics_process.html)     | [condition_non_historical](../Domains/Temporal/temporal_containers.html#temporal-condition-non-historical)               |            |
| [Distribution](../Statistical_Methods/Distribution.html)         | [Trace](../Norm_Methods/trace.html)         |                                                               | [Temporal statistics process](../Domains/Temporal/temporal_statistics_process.html)                  |            |
| [Norm methods](../Norm_Methods/Overview.html)         | [Index](../Norm_Methods/index_based.html)         |                                                               |                                                                              |            |
|                                       | [Component](../Norm_Methods/component_based.html) |                                                               |                                                                              |            |

If you prefer to use statistics of a simulation using StatisticsApplication, there is `spatial_statistics_process` for spatial statistics calculations and `temporal_statistics_process` for temporal statistics. These processes can be included via JSON settings under `auxiliary_processes`.