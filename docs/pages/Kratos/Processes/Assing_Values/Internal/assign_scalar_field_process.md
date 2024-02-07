---
title: Assign Scalar Field
keywords: process core
tags: [assign scalar field to entities process]
sidebar: kratos_core_processes
summary: 
---

# Assign Scalar Field

## Description

This process assigns a value depending on a function to a variable in all of the entities in a given `ModelPart`.

This process has several entry points:
__Python__:
- `AssignScalarInputToNodesProcess`
- `AssignScalarInputHistoricalToNodesProcess`
- `AssignScalarInputToConditionsProcess`
- `AssignScalarInputToElementsProcess`

__CPP__:
- `AssignScalarInputToEntitiesProcess<Entity>`

Depending on the type of the target variable, the behaviour is the following:
- `Variable<double>`: It is evaluated in the center of the entities
- `array_1d_component_type`: The same as Variable evaluated in the center of the element
- `Variable<Vector>`: The vector has to have a size equal to the number of nodes, and its values are computed per each entry of the vector using the coordinates of the nodes which occupies the same position in the geometry

This process is executed on the follwing hooks:
- `ExecuteInitializeSolutionStep`

## Parameters & Defaults

```json
{
    "model_part_name" :"model_part_name",
    "mesh_id"         : 0,
    "variable_name"   : "variable_name",
    "interval"        : [0.0, 1e30],
    "value"           : "please give an expression in terms of the variable x, y, z, t",
    "local_axes"      : {}
}
```

`model_part_name`: Name of the modelpart in wich the field variable will be applied
`mesh_id`: Id of the internal mesh to which the process will be applied. Default `0`.
`variable_name`: Name of the variable in which the field value will be applied.
`interval`: Interval of time in which the process will be applied.
`value`: Expression that will be evaluated for every 
`local_axes`: 
