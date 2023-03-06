---
title: Temporal Methods
keywords: statistics, temporal, methods
tags: [temporal_methods.md]
sidebar: statistics_application
summary: 
---
All the temporal methods are available through `TemporalMethods` submodule under `StatisticsApplication`. In the case of temporal methods, seperate objects needs to be created for each input variable for different type of variables as shown in following example for nodal non historical data container value methods.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_double_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Double(model_part, "", Kratos.PRESSURE, 0, Kratos.DENSITY)
sum_array_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_vector_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Vector(model_part, "", Kratos.LOAD_MESHES, 0, Kratos.MATERIAL_PARAMETERS)
sum_matrix_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Matrix(model_part, "", Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, 0, Kratos.CAUCHY_STRESS_TENSOR)
```

In order to remove the hassle of using different names (such as `Double`, `Array`...), one can get the advantage of common methods available in each of the containers. Following example will create the same objects as in the previous example in one go using simplified sum method.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
params = Kratos.Parameters(r"""
{
    "input_variables" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
    "output_variables" : ["DENSITY", "DISPLACEMENT", "MATERIAL_PARAMETERS", "CAUCHY_STRESS_TENSOR"]
}""")
sum_methods = KratosStats.TemporalMethods.NonHistorical.Nodes.Sum(model_part, "none", 0, params)
```

For **Sum**, **Mean**, **RootMeanSquare** methods, above `Kratos.Parameters` object keys will remain the same. But for **Variance** method following `json` parameters are used.

```json
{
    "input_variables"           : [],
    "output_mean_variables"     : [],
    "output_variance_variables" : []
}
```

When using these simplified, one needs to take care of compatibility of input and output variables depending on the method and norm type, otherwise runtime errors will be thrown.