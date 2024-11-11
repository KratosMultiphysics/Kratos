---
title: Assign Vector By Direction
keywords: process core
tags: [assign Vector variable by direction]
sidebar: kratos_core_processes
summary: 
---

# Assign Vector Variable By Direction

## Description

This process assigns a value to a variable in all of the entities in a given `ModelPart`.

This process has a generic versions for vectors:
- `AssignVectorByDirectionToEntityProcess`

Different specialized entry points for directional vectors
- `AssignVectorByDirectionProcess`
- `AssignVectorByDirectionToNodeProcess`
- `AssignVectorByDirectionToElementProcess`
- `AssignVectorByDirectionToConditionProcess`

This process will call the corresponding `AssignScalarVariableToEntityProcess` for each component of the vector.

## Execution

This process is executed in the following hooks:

#### `ExecuteInitializeSolutionStep`

Assigns the value. If `constrained` is `true` fixes the values after assigning them.

#### `ExecuteFinalizeSolutionStep`

If `constrained` is `true` unfixes the values.

## Parameters & Defaults

```json
{
    "model_part_name" : "please_specify_model_part_name",
    "variable_name"   : "SPECIFY_VARIABLE_NAME",
    "interval"        : [0.0, 1e30],
    "modulus"         : 1.0,
    "direction"       : [1.0, 0.0, 0.0],
    "constrained"     : [true,true,true],
    "local_axes"      : {},
    "entities"        : []
}
```

##### `model_part_name` 
Name of the modelpart in wich the process will be applied.

##### `variable_name`
Name of the variable in which the process will be applied.

##### `interval`
Specifies the interval of time in which the process will be applied.

##### `modulus`
Scalar value that will define the modulus of the direction vector. 

##### `direction`
Vector that will define the of direction. It can be set to `"Automatic"` 

##### `constrained`
Specifies if the variable has a constrain `true` or not `false` over each dimension.

Setting to `true` will cause the values to be fixed during `ExecuteInitializeSolutionStep` and unfixed during `ExecuteFinalizeSolutionStep`.

##### `local_axes`
Specifies the axes in which the value will be applied

##### `entities`
This option is only available for the `AssignVectorVariableToEntitiesProcess` entry Point.

List of entities into which the value will be applies. Accepts: `nodes`, `elements`, `conditions`.
