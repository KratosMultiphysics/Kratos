---
title: Python Modules For Processes
keywords: json
tags: [PythonModulesForProcesses.md]
sidebar: structural_mechanics_application
summary:
---
# Python Modules for Processes

## AssignVectorVariableProcess:
- This process assigns a given value (vector) to the nodes belonging to a certain submodelpart.
- Following are the default settings for the parameters of this process inside the ProjectParameters.json

```json
"Parameters"    :{
    "mesh_id"              : 0,
    "model_part_name"      : "please_specify_model_part_name",
    "variable_name"        : "SPECIFY_VARIABLE_NAME",
    "interval"             : [0.0, 1e30],
    "value"                : [0.0, 0.0, 0.0],
    "constrained"          : [true,true,true],
    "local_axes"           : {}
}
```
- example of admissible values for "value" : [10.0, "3*t", "x+y"]

## AssignVectorByDirectionToConditionProcess:
- This process sets a variable with a certain scalar value in a given direction, for all the 'conditions' belonging to a submodelpart.
- Following are the default settings for the parameters of this process inside the ProjectParameters.json

```json
"Parameters"    :{
    "mesh_id"              : 0,
    "model_part_name"      : "please_specify_model_part_name",
    "variable_name"        : "SPECIFY_VARIABLE_NAME",
    "interval"             : [0.0, 1e30],
    "modulus"              : 0.0,
    "direction"            : [1.0, 0.0, 0.0],
    "local_axes"           : {},
    "entities"             : ["conditions"]
}
```
- The value can be a double or a string.

## AssignVectorByDirectionProcess:
- This process sets a variable with a certain scalar value in a given direction, for all the 'nodes' belonging to a submodelpart.
- Following are the default settings for the parameters of this process inside the ProjectParameters.json

```json
"Parameters"    :{
    "mesh_id"              : 0,
    "model_part_name"      : "please_specify_model_part_name",
    "variable_name"        : "SPECIFY_VARIABLE_NAME",
    "interval"             : [0.0, 1e30],
    "modulus"              : 0.0,
    "constrained"          : true,
    "direction"            : [1.0, 0.0, 0.0],
    "local_axes"           : {}
}
```

## AssignScalarVariableToConditionsProcess:
- This process sets a variable a certain scalar value, for all the conditions belonging to a submodelpart. As this is used for scalar variables, it does not need a direction.
- Following are the default settings for the parameters of this process inside the ProjectParameters.json

```json
"Parameters"    :{
    "mesh_id"         : 0,
    "model_part_name" : "please_specify_model_part_name",
    "variable_name"   : "SPECIFY_VARIABLE_NAME",
    "interval"        : [0.0, 1e30],
    "value"           : 0.0,
    "local_axes"      : {},
    "entities"        : ["conditions"]
}
```

## process_factory:
- If a python_module is not available for the process you want to use, then use process_factory; eg. for ApplyConstantVectorValueProcess.
- Example:

```json
"loads_process_list"       : [{
    "python_module"   : "process_factory",
    "kratos_module" : "KratosMultiphysics",
    "process_name"          : "ApplyConstantVectorValueProcess",
    "Parameters"            : {
        "mesh_id"         : 0,
        "model_part_name" : "Structure.PointLoad3D_load",
        "variable_name"   : "POINT_LOAD",
        "modulus"         :  1,
        "direction"       : [0.0,0.0,1.0]
    }
}]
```