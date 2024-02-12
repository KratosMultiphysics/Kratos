---
title: Check Scalar
keywords: process core
tags: [Check Scalar Base]
sidebar: kratos_core_processes
summary: 
---

# Check Scalar

## Description

This process is the base class to check analytically from a function the solution (scalar) in a certain entity belonging a certain submodelpart.

This process has the folling entry points:

- `CheckScalarBaseProcess`: Base process

- `CheckScalarFromProcessInfoProcess`: Checks analytically from a function the solution (scalar) in the process info belonging a certain submodelpart

- `CheckScalarToNodesProcess`: Checks analytically from a function the solution (scalar) in a set of nodes belonging a certain submodelpart

## Parameters & Defaults

```json
{
    "model_part_name" : "please_specify_model_part_name",
    "mesh_id"         : 0,
    "variable_name"   : "SPECIFY_VARIABLE_NAME",
    "interval"        : [0.0, 1e30],
    "value"           : 0.0,
    "tolerance_rank"  : 3
}
```

##### `model_part_name` 
Name of the modelpart to check.

##### `mesh_id` 
Id of the msh to be check.

##### `variable_name` 
Name of the variable to check.

##### `interval` 
Interval of time in which the check will be executed.

##### `value` 
Reference value to check. It can be a string with a mathematical expression containin the following variables: `x`, `y`, `z`, `t`.

##### `tolerance_rank` 
Tolerance of the value checked.
