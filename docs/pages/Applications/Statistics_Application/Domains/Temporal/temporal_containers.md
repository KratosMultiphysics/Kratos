---
title: Temporal Containers
keywords: statistics, temporal, containers
tags: [temporal_containers.md]
sidebar: statistics_application
summary: 
---
In the temporal domain, there are five different containers available for user to calculate temporal statistics. They are,

- [Temporal nodal historical historical](#temporal-nodal-historical-historical)
- [Temporal nodal historical non historical](#temporal-nodal-historical-non-historical)
- [Temporal nodal non historical](#temporal-nodal-non-historical)
- [Temporal condition non historical](#temporal-condition-non-historical)
- [Temporal element non historical](#temporal-element-non-historical)

#### Temporal nodal historical historical

In this container, statistics are calculated on user specified variables in the nodal historical container and, output is also written to the nodal historical container. Example is given below.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.TemporalMethods.Historical.HistoricalOutput.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_norm = KratosStats.TemporalMethods.Historical.HistoricalOutput.NormMethods.Sum.Array(model_part, "magnitude", Kratos.VELOCITY, 0, Kratos.DENSITY)
```

#### Temporal nodal historical non historical

In this container, statistics are calculated on user specified variables in the nodal historical container and, output is also written to the nodal non historical container. Example is given below.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.TemporalMethods.Historical.NonHistoricalOutput.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_norm = KratosStats.TemporalMethods.Historical.NonHistoricalOutput.NormMethods.Sum.Array(model_part, "magnitude", Kratos.VELOCITY, 0, Kratos.DENSITY)
```

#### Temporal nodal non historical

In this container, statistics are calculated on user specified variables in the nodal non historical container and, output is also written to the nodal non historical container. Example is given below.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_norm = KratosStats.TemporalMethods.NonHistorical.Nodes.NormMethods.Sum.Array(model_part, "magnitude", Kratos.VELOCITY, 0, Kratos.DENSITY)
```

#### Temporal condition non historical

In this container, statistics are calculated on user specified variables in the condition non historical container and, output is also written to the condition non historical container. Example is given below.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.TemporalMethods.NonHistorical.Conditions.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_norm = KratosStats.TemporalMethods.NonHistorical.Conditions.NormMethods.Sum.Array(model_part, "magnitude", Kratos.VELOCITY, 0, Kratos.DENSITY)
```

#### Temporal element non historical

In this container, statistics are calculated on user specified variables in the element non historical container and, output is also written to the element non historical container. Example is given below.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.TemporalMethods.NonHistorical.Elements.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_norm = KratosStats.TemporalMethods.NonHistorical.Elements.NormMethods.Sum.Array(model_part, "magnitude", Kratos.VELOCITY, 0, Kratos.DENSITY)
```