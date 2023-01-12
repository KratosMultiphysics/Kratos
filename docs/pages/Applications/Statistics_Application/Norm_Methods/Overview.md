---
title: Overview
keywords: norm, statistics
tags: [Overview.md]
sidebar: statistics_application
summary: 
---

## Norms Summary

Few different norms are predefined in this application. Following table summarize supported variable types for each norm.

|                                   | Double | Array 3D | Vector | Matrix |
|-----------------------------------|--------|----------|--------|--------|
| [Value](value.html)               | [x]    |          |        |        |
| [Magnitude](magnitude.html)       | [x]    | [x]      | [x]    | [x]    |
| [Euclidean](euclidean.html)       |        | [x]      | [x]    |        |
| [Infinity](infinity.html)         |        | [x]      | [x]    | [x]    |
| [P-Norm](p_norm.html)             |        | [x]      | [x]    | [x]    |
| [Lpq-Norm](lpq_norm.html)         |        |          |        | [x]    |
| [Frobenius](frobenius.html)       |        |          |        | [x]    |
| [Trace](trace.html)               |        |          |        | [x]    |
| [Index](index_based.html)         |        |          | [x]    | [x]    |
| [Component](component_based.html) |        | [x]      |        |        |

## Supported methods summary

All of the value methods mentioned in ```Statistical Methods``` supports norm version of it, in these methods, higher dimensional values are transformed in to scalar values by a use specified norm, and then statistics are calculated based on the chosen method. Supported norm types may differ based on the variable data type being used. There are few methods, which only supports norm methods. Following table summarize availability of value and norm methods.

| Statistical Methods                   | Value Method | Norm Method |
|---------------------------------------|--------------|-------------|
| [Sum](../Statistical_Methods/Sum.html)                           | [x]          | [x]         |
| [Mean](../Statistical_Methods/Mean.html)                         | [x]          | [x]         |
| [Root mean square](../Statistical_Methods/Root_Mean_Square.html) | [x]          | [x]         |
| [Variance](../Statistical_Methods/Variance.html)                 | [x]          | [x]         |
| [Min](../Statistical_Methods/Min.html)                           |              | [x]         |
| [Max](../Statistical_Methods/Max.html)                           |              | [x]         |
| [Median](../Statistical_Methods/Median.html)                     |              | [x]         |
| [Distribution](../Statistical_Methods/Distribution.html)         |              | [x]         |

## Example

Following example shows variance method, under **value** category and **norm** category for nodal non historical velocity. Norm method uses `"magnitude"` as the norm to reduce `VELOCITY` to a scalar value. `value_mean` and `value_variance` will be of `Array 3D` type, whereas `norm_mean` and `norm_variance` will be of `Double` type.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
value_mean, value_variance = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Variance(model_part, Kratos.VELOCITY)
norm_mean, norm_variance = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Variance(model_part, Kratos.VELOCITY, "magnitude")
```