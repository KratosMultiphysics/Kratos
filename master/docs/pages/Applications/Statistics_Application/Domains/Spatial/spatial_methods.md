---
title: Spatial Methods
keywords: statistics, spatial, methods
tags: [spatial_methods.md]
sidebar: statistics_application
summary: 
---
All the spatial statistical methods can be found in the `SpatialMethods` submodule. `ValueMethods` container will contain all available spatial value methods, and `NormMethods` container will contain all available norm methods.

Following example illustrates all the available value methods executed on nodal non historical variable velocity.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Sum(model_part, Kratos.VELOCITY)
mean_value = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Mean(model_part, Kratos.VELOCITY)
rms_value = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.RootMeanSquare(model_part, Kratos.VELOCITY)
mean_value, variance_value = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Variance(model_part, Kratos.VELOCITY)
```

Following example illustrates all the available norm methods executed on nodal non historical variable velocity's magnitude.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_norm = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Sum(model_part, Kratos.VELOCITY, "magnitude")
mean_norm = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Mean(model_part, Kratos.VELOCITY, "magnitude")
rms_norm = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.RootMeanSquare(model_part, Kratos.VELOCITY, "magnitude")
mean_norm, variance_norm = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Variance(model_part, Kratos.VELOCITY, "magnitude")
min_norm, min_norm_id = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Min(model_part, Kratos.VELOCITY, "magnitude")
max_norm, max_norm_id = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Max(model_part, Kratos.VELOCITY, "magnitude")
min_norm, max_norm, group_upper_norms, group_histogram, group_percentage_distribution = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Distribution(model_part, Kratos.VELOCITY, "magnitude")
```