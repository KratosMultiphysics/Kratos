---
title: Point, Line and Surface Loads
keywords: mpm load
tags: [mpm load]
sidebar: mpm_application
summary: 
---

The [`AssignVectorByDirectionToConditionProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/assign_vector_by_direction_to_condition_process.py) is used in the
`MPMApplication` to impose point, line and surface loads to the conditions of a background grid submodel part.

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics",
    "python_module" : "assign_vector_by_direction_to_condition_process",
    "Parameters"    : {
        "model_part_name"      : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "variable_name"        : "SPECIFY_VARIABLE_NAME",
        "interval"             : [0.0, 1e30],
        "modulus"              : 0.0,
        "direction"            : [1.0, 0.0, 0.0],
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart.

##### `variable_name`
String identifying the variable to which the vector value has to be assigned.
Admissible values are:
* `POINT_LOAD`
* `LINE_LOAD`
* `SURFACE_LOAD`

##### `interval`
Time interval in which the process applies.

##### `modulus`
Modulus of the vector to be assigned to the variable specified by `variable_name`.
The value can be either a number or a string representing a function depending on
space and/or time.

##### `direction`
Direction of the vector to be assigned to the variable specified by `variable_name`. Each component of the list can be either a number or a string representing a function depending on space and/or time.

## Source Code

[<i class="fa fa-github"></i> `kratos/python_scripts/assign_vector_by_direction_to_condition_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/assign_vector_by_direction_to_condition_process.py)
