---
title: Variance
keywords: 
tags: [Variance.md]
sidebar: statistics_application
summary: 
---
## Definition

In the case of spatial domain, it calculates variance of a given variable for a given container and returns mean and variance as shown in following equations. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Results will have the same type as the type of the variable specified by the user. (If it has higher dimension than a scalar, mean of each dimension will be calculated seperately resulting with a mean having same dimension as the input dimension)

<p align="center">$$\underline{\bar{x}} = \frac{1}{N}\sum_{i=1}^N{\underline{x}_i} $$</p>
<p align="center">$$\underline{v} = \frac{1}{N}\sum_{i=1}^N{\left(\underline{x}_i - \underline{\bar{x}} \right )^2}$$</p>

In the case of temporal domain, **Variance** method is the time integrated quantity's variance for a specific variable. Mean and variance will be stored each element under user specified variables and a user specified container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Results will have the same type as the type of the variable specified by the user preserving the dimensionality as in the spatial case.

<p align="center">$$ \underline{\bar{x}} = \frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k\Delta t_k}$$</p>
<p align="center">$$Var\left(\underline{x} \right ) = \frac{1}{T_{total}}\sum_{k=1}^{P}{\left(\underline{x}_k - \underline{\bar{x}} \right )^2\Delta t_k}$$</p>
<p align="center">$$ \quad where $$</p>
<p align="center">$$T_{total} = T_{end} - T_{initial} \quad and \quad \Delta t_k = T_{k} - T_{k-1} \quad \forall T_k \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace$$</p>

## Examples

### Spatial

Following is an example of variance calculation of non historical `VELOCITY` over the whole model part's nodes

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
mean, variance = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Variance(model_part, Kratos.VELOCITY)
```

### Temporal

Following is an example of root mean square calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable `VECTOR_3D_MEAN` will store mean and `VECTOR_3D_VARIANCE` will store variance in same container for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
variance_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Variance.Array(model_part, "", Kratos.VELOCITY, 0, KratosStats.VECTOR_3D_MEAN, KratosStats.VECTOR_3D_VARIANCE)
integration_starting_time = 2.0
variance_method.InitializeStatisticsMethod(integration_starting_time)
for t in range(3, 6):
    variance_method.CalculateStatistics()
```