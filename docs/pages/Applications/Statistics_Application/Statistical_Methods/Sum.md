---
title: Sum
keywords: 
tags: [Sum.md]
sidebar: statistics_application
summary: 
---
## Definition

In the case of spatial domain, it adds up all the variable values for a given container and returns summed up value as shown in following equation. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user

<p align="center">$$ \underline{r} = \sum_{i=1}^N{\underline{x}_i}$$</p>

In the case of temporal domain, **Sum** methods is the time integrated quantity for a specific variable. It will be stored each element under user specified variable and a user specified container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user.

<p align="center">$$ \underline{r} = \sum_{k=1}^{P}{\underline{x}_k\Delta t_k} \quad where \quad \Delta t_k = T_{k} - T_{k-1} \quad \forall T_k \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace$$</p>

## Examples

### Spatial

Following is an example of summation of non historical `VELOCITY` over the whole model part's nodes

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Sum(model_part, Kratos.VELOCITY)
```

### Temporal

Following is an example of integration calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable is same containers `DISPLACEMENT` where integrated value will be stored for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
integration_starting_time = 2.0
sum_method.InitializeStatisticsMethod(integration_starting_time)
for t in range(3, 6):
    sum_method.CalculateStatistics()
```