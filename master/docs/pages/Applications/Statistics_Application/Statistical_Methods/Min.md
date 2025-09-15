---
title: Min
keywords: 
tags: [Min.md]
sidebar: statistics_application
summary: 
---
## Definition

In the case of spatial domain, it returns minimum of a given variable's norm for a given container, and the corresponding items id. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Results will be double and integer, irrespective of the input type, since higher dimensional variable types will be reduced to scalars by the use of norms.

<p align="center">$$v = \min_{\underline{x}_i \in \mathbf{T}} |\underline{x}_i|$$</p>

In the case of temporal domain, **Min** method returns minimum value in the temporal domain, and the time minimum is found. Minimum and its occurring time will be stored each element under user specified variables and a user specified container. x<sub>k</sub> is the k<sup>th</sup> time step element's variable value of the corresponding container. Results will have the same type as the type of the variable specified by the user preserving the dimensionality as in the spatial case.

<p align="center">$$v = \min_{\underline{x}_k \in \mathbf{T}} |\underline{x}_k|$$</p>

## Examples

### Spatial

Following is an example of min method of non historical `VELOCITY`'s magnitude over the whole model part's nodes. It returns a tuple, first argument being the minimum, and the second argument being the id of the node where the minimum is found.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
min_value, min_id = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Min(model_part, Kratos.VELOCITY, "magnitude")
```

### Temporal

Following is an example of min method in non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable `VECTOR_3D_NORM` will store minimum and `TIME` will store the time minimum occured for each node. The `0` represents echo level for this method object. "magnitude" indicates that magnitude norm is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
min_method = KratosStats.TemporalMethods.NonHistorical.Nodes.NormMethods.Min.Array(model_part, "magnitude", Kratos.VELOCITY, 0, KratosStats.VECTOR_3D_NORM, Kratos.TIME)
integration_starting_time = 2.0
min_method.InitializeStatisticsMethod(integration_starting_time)
for t in range(3, 6):
    min_method.CalculateStatistics()
```