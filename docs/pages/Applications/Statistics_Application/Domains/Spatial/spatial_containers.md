---
title: Spatial Containers
keywords: statistics, spatial, containers
tags: [spatial_containers.md]
sidebar: statistics_application
summary: 
---
Four different types of containers are supported for spatial methods.

- [Spatial nodal historical](#spatial-nodal-historical)
- [Spatial nodal non historical](#spatial-nodal-non-historical)
- [Spatial condition non historical](#spatial-condition-non-historical)
- [Spatial element non historical](#spatial-element-non-historical)

#### Spatial nodal historical

Nodal historical containers methods can be accessed through `Historical` submodule in `SpatialMethods`. This calculates chosen statistics of the variable in the nodal historical data. Following example illustrates use of value and norm methods.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
sum_value = KratosStats.SpatialMethods.Historical.ValueMethods.Sum(model_part, Kratos.VELOCITY)
sum_norm = KratosStats.SpatialMethods.Historical.NormMethods.Sum(model_part, Kratos.VELOCITY, "magnitude")
```

#### Spatial nodal non historical

Nodal non historical containers methods can be accessed through `NonHistorical.Nodes` submodule in `SpatialMethods`. This calculates chosen statistics of the variable in the nodal non historical data. Following example illustrates use of value and norm methods.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Sum(model_part, Kratos.VELOCITY)
sum_norm = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Sum(model_part, Kratos.VELOCITY, "magnitude")
```

#### Spatial condition non historical

Condition non historical containers methods can be accessed through `NonHistorical.Conditions` submodule in `SpatialMethods`. This calculates chosen statistics of the variable in the condition non historical data. Following example illustrates use of value and norm methods.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.SpatialMethods.NonHistorical.Conditions.ValueMethods.Sum(model_part, Kratos.VELOCITY)
sum_norm = KratosStats.SpatialMethods.NonHistorical.Conditions.NormMethods.Sum(model_part, Kratos.VELOCITY, "magnitude")
```

#### Spatial element non historical

Element non historical containers methods can be accessed through `NonHistorical.Elements` submodule in `SpatialMethods`. This calculates chosen statistics of the variable in the element non historical data. Following example illustrates use of value and norm methods.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.SpatialMethods.NonHistorical.Elements.ValueMethods.Sum(model_part, Kratos.VELOCITY)
sum_norm = KratosStats.SpatialMethods.NonHistorical.Elements.NormMethods.Sum(model_part, Kratos.VELOCITY, "magnitude")
```