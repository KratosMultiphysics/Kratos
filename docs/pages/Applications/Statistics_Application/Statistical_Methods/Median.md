---
title: Median
keywords: 
tags: [Median.md]
sidebar: statistics_application
summary: 
---
## Definition

Median returns the median value in spatial domain of a given variable's norm for a given container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Results will be double value type, irrespective of the input type, since higher dimensional variable types will be reduced to scalars by the use of norms.

## Example

Following is an example of median method of non historical `VELOCITY`'s magnitude over the whole model part's nodes. It returns a double which is the median.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
median_value = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Median(model_part, Kratos.VELOCITY, "magnitude")
```