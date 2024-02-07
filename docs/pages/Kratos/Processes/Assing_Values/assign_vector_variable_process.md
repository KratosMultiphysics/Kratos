---
title: Assign Scalar Variable
keywords: process core
tags: [assign scalar variable to entities process]
sidebar: kratos_core_processes
summary: 
---

# Assign Scalar Variable

## Description

This process assigns a value to a variable in all of the entities in a given `ModelPart`.

This process has a generic version:
- `AssignScalarVariableToEntitiesProcess`

Different specialized entry points for scalar
- `AssignScalarVariableToProcess`
- `AssignScalarVariableToNodesProcess`
- `AssignScalarVariableToElementsProcess`
- `AssignScalarVariableToConditionsProcess`

## Execution

This process is executed in the following hooks:

#### `ExecuteInitializeSolutionStep`

Assigns the value.


## Parameters & Defaults

```json
{
    "model_part_name" :"model_part_name",
    "mesh_id"         : 0,
    "variable_name"   : "variable_name",
    "interval"        : [0.0, 1e30],
    "value"           : [1.0, 1.0, 1.0],
    "local_axes"      : {},
    "entities"        : []
}
```

##### `model_part_name` 
Name of the modelpart in wich the process will be applied.

##### `mesh_id`
Id of the internal mesh to which the process will be applied. Default `0`.

##### `variable_name`
Name of the variable in which the process will be applied.

##### `interval`
Specifies the interval of time in which the process will be applied.

##### `value`
Vector value that will be applied. 

It can be expressed as an array of `[double]` or `[string]` containing a mathematical expression with `x`,`y`,`z`,`t` as variables for each position of the array.

##### `local_axes`
Specifies the axes in which the value will be applied

##### `entities`
This option is only available for the `AssignScalarVariableToEntitiesProcess` entry Point.

List of entities into which the value will be applies. Accepts: `nodes`, `elements`, `conditions`.