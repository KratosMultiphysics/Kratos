---
title: Importing the Kratos Library for Scalar Assignments
keywords: process core, Kratos library
tags: [assign scalar variable to entities]
sidebar: kratos_core_processes
summary: This document provides an overview of a process in the Kratos Multiphysics library that assigns a scalar variable to various entities within a model part.
---

# Assign Scalar Variable To Entities Process

## Description

This process assigns a predefined scalar value to specified entities within a given `ModelPart` in Kratos Multiphysics simulations.

This process has specific versions depending on the entities being targeted:
- `AssignScalarVariableToNodesProcess`
- `AssignScalarVariableToElementsProcess`
- `AssignScalarVariableToConditionsProcess`
- `AssignScalarVariableToMasterSlaveConstraintsProcess`

The scalar value assigned can be either a constant number or a function defined as a string, which is evaluated at runtime.

## Execution

This process is executed during the following hooks:

#### `ExecuteInitializeSolutionStep`

During this step, the scalar value is assigned to the entities. If any conditions are specified, they are also applied.

#### `ExecuteFinalizeSolutionStep`

Finalizes any changes or updates made during the initialization step, ensuring the state is consistent moving forward.

## Parameters & Defaults

```json
{
    "model_part_name" : "please_specify_model_part_name",
    "variable_name"   : "SPECIFY_VARIABLE_NAME",
    "interval"        : [0.0, 1e30],
    "value"           : 0.0,
    "local_axes"      : {},
    "historical"      : false,
    "entities"        : []
}
```

##### `model_part_name` 
Defines the model part where the scalar variable will be applied.

##### `variable_name`
Specifies the name of the scalar variable to be assigned.

##### `interval`
Indicates the time interval over which the scalar assignment is active.

##### `value`
Determines the scalar value or function to be assigned to the variable.

##### `local_axes`
Defines the local axes configuration for the assignment if needed.

##### `historical`
Specifies whether the assignment should be stored historically.

##### `entities`
Lists the types of entities to which the scalar value will be assigned. Options include `nodes`, `elements`, `conditions`, and `constraints`. This is essential as it dictates the scope of the process.