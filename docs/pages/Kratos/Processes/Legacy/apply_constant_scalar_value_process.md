---
title: Apply Constant Scalarvalue
keywords: process core
tags: [apply constant scalarvalue process]
sidebar: kratos_core_processes
summary: 
---

# Apply Constant Scalarvalue

## Warning

This is a Legacy process and should not be used directly. Please see [Assign Scalar Variable To Entities](/Assign_Values/assign_scalar_variable_to_entities_process.md)

## Description

This process is intended to be used for setting scalar values (`bool`, `int`, `double`) to a variable in a model part.

This process is executed on the follwing hooks:
- `ExecuteInitialize`

## Parameters & Defaults

```json
{
    "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
    "mesh_id": 0,
    "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
    "is_fixed": false,
    "value" : 1.0
}
```

`model_part_name`: Name of the modelpart in wich the scalar variable will be applied
`mesh_id`: Id of the internal mesh to which the process will be applied. Default `0`.
`variable_name`: Name of the variable in which the scalar value will be applied.
`is_fixed`: Sets the value as fixed (`true`) or not fixed (`false`).
`value`: Scalar value that will be applied.
