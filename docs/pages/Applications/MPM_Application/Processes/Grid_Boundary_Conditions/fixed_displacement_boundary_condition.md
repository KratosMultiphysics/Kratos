---
title: Fixed Displacement
keywords: mpm fixed displacement
tags: [mpm fixed displacement]
sidebar: mpm_application
summary: 
---

A zero displacement Dirichlet boundary condition can be imposed on the nodes of a given background grid submodel part by means of the python process [`AssignVectorVariableProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/assign_vector_variable_process.py) implemented in the Kratos Core.

## Parameters & Defaults

The parameters to be used are the following:

```json
{
    "kratos_module" : "KratosMultiphysics",
    "python_module" : "assign_vector_variable_process",
    "Parameters"    : {
        "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "variable_name"   : "DISPLACEMENT",
        "interval"        : [0.0, 1e30],
        "constrained"     : [true, true, true],
        "value"           : [0.0, 0.0, 0.0]
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart.

##### `interval`
Time interval in which the process applies.

## Process File

[<i class="fa fa-github"></i> `kratos/python_scripts/assign_vector_variable_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/assign_vector_variable_process.py)
