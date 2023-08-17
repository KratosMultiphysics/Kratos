---
title: Mean
keywords: 
tags: [Mean.md]
sidebar: statistics_application
summary: 
---
## Definition

In the case of spatial domain, it calculates mean of a given variable for a given container and returns it as shown in following equation. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user. (If it has higher dimension than a scalar, mean of each dimension will be calculated seperately resulting with a mean having same dimension as the input dimension)

<p align="center">$$ \underline{r} = \frac{1}{N}\sum_{i=1}^{N}\underline{x}_i$$</p>


In the case of temporal domain, **Mean** methods is the time integrated quantity's mean for a specific variable. It will be stored each element under user specified variable and a user specified container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user preserving the dimensionality as in the spatial case.

<p align="center">$$\underline{\bar{x}} = \frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k\Delta t_k} \quad where \quad T_{total} = T_{end} - T_{initial} \quad and \quad \Delta t_k = T_{k} - T_{k-1} \quad \forall T_k \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace $$</p>

## Examples

### Spatial

Following is an example of mean calculation of non historical `VELOCITY` over the whole model part's nodes

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
mean = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Mean(model_part, Kratos.VELOCITY)
```

### Temporal

Following is an example of mean calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable is same containers `VECTOR_3D_MEAN` where mean will be stored for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
mean_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Mean.Array(model_part, "", Kratos.VELOCITY, 0, KratosStats.VECTOR_3D_MEAN)
integration_starting_time = 2.0
mean_method.InitializeStatisticsMethod(integration_starting_time)
for t in range(3, 6):
    mean_method.CalculateStatistics()
```
