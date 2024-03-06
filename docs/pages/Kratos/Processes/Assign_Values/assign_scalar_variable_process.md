---
title: Assign Scalar Variable
keywords: process core
tags: [assign scalar variable to entities process]
sidebar: kratos_core_processes
summary: 
---

# Assign Scalar Variable

## Description

This process assigns a scalar value to a variable in all of the entities in a given `ModelPart`.

For the vector version please see [Assign Vector Variable](/Assign_Values/assign_vector_variable_process)

This process has a generic version:
- `AssignScalarVariableToEntitiesProcess`

Different specialized entry points for scalars
- `AssignScalarVariableToProcess`
- `AssignScalarVariableToNodesProcess`
- `AssignScalarVariableToElementsProcess`
- `AssignScalarVariableToConditionsProcess`
- `AssignScalarVariableToConstraintsProcess`

Different specialized entry points for fields
- `AssignScalarFieldToNodesProcess`
- `AssignScalarFieldElementsProcess`
- `AssignScalarFieldConditionsProcess`

## Execution

This process is executed in the following hooks:

#### `ExecuteInitializeSolutionStep`

Assigns the value. If `constrained` is `true` fixes the values after assigning them.

#### `ExecuteFinalizeSolutionStep`

If `constrained` is `true` unfixes the values.

## Parameters & Defaults

```json
{
    "model_part_name" :"model_part_name",
    "mesh_id"         : 0,
    "variable_name"   : "variable_name",
    "interval"        : [0.0, 1e30],
    "value"           : 1.0,
    "constrained"     : true,
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
Scalar value that will be applied. 

It can be expressed as a `bool`, `int`, `double` or `string` containing a mathematical expression with `x`,`y`,`z`,`t` as variables.

##### `constrained`
Specifies if the variable has a constrain `true` or not `false`.

Setting to `true` will cause the values to be fixed during `ExecuteInitializeSolutionStep` and unfixed during `ExecuteFinalizeSolutionStep`.

##### `local_axes`
Specifies the axes in which the value will be applied

##### `entities`
This option is only available for the `AssignScalarVariableToEntitiesProcess` entry Point.

List of entities into which the value will be applies. Accepts: `nodes`, `elements`, `conditions`, `constraints`.