---
title: Check Vector Components To Nodes
keywords: process core
tags: [Check Vector Components To Nodes]
sidebar: kratos_core_processes
summary: 
---

# Check Vector Components To Nodes

## Description

This process checks analytically from a function the solution (vector) in a set of nodes belonging a certain submodelpart

For every vector component this process will call [`checl_scalar_process`]()

## Parameters & Defaults

```json
{
    "model_part_name" : "please_specify_model_part_name",
    "mesh_id"         : 0,
    "variable_name"   : "SPECIFY_VARIABLE_NAME",
    "interval"        : [0.0, 1e30],
    "value"           : [10.0, "3*t", "x+y"],
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
Reference value array to check. Each possition can be double or a string with a mathematical expression containin the following variables: `x`, `y`, `z`, `t`.

##### `tolerance_rank` 
Tolerance of the value checked.
