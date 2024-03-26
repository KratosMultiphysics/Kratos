---
title: Python Modules For Processes
keywords: json
tags: [PythonModulesForProcesses.md]
sidebar: structural_mechanics_application
summary:
---
# Python Modules for Processes

## AssignVectorVariableProcess:
- This process assigns a given value (vector) to the nodes belonging to a certain modelpart.
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
- Sources:
    * [kratos/python_scripts/assign_vector_variable_process.py](https://github.com/KratosMultiphysics/Kratos/tree/master/kratos/python_scripts/assign_vector_variable_process.py)
    * [Doxygen](https://kratos-docs.onrender.com/d7/dd4/classassign__vector__variable__process_1_1_assign_vector_variable_process.html)

## AssignVectorByDirectionToConditionProcess:
- This process sets a variable with a certain scalar value in a given direction, for all the 'conditions' belonging to a modelpart.
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
- Sources:
    * [kratos/python_scripts/assign_vector_by_direction_to_condition_process.py](https://github.com/KratosMultiphysics/Kratos/tree/master/kratos/python_scripts/assign_vector_by_direction_to_condition_process.py)
    * [Doxygen](https://kratos-docs.onrender.com/d7/d98/classassign__vector__by__direction__to__condition__process_1_1_assign_vector_by_direction_to_condition_process.html)

## AssignVectorByDirectionProcess:
- This process sets a variable with a certain scalar value in a given direction, for all the 'nodes' belonging to a modelpart.
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
- Sources:
    * [kratos/python_scripts/assign_vector_by_direction_process.py](https://github.com/KratosMultiphysics/Kratos/tree/master/kratos/python_scripts/assign_vector_by_direction_process.py)
    * [Doxygen](https://kratos-docs.onrender.com/da/d16/classassign__vector__by__direction__process_1_1_assign_vector_by_direction_process.html)

## AssignScalarVariableToConditionsProcess:
- This process sets a variable with a certain scalar value, for all the conditions belonging to a modelpart. As this process is used for scalar variables, it does not need a direction.
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
- Sources:
    * [kratos/python_scripts/assign_scalar_variable_to_conditions_process.py](https://github.com/KratosMultiphysics/Kratos/tree/master/kratos/python_scripts/assign_scalar_variable_to_conditions_process.py)
    * [Doxygen](https://kratos-docs.onrender.com/d3/da3/classassign__scalar__variable__to__conditions__process_1_1_assign_scalar_variable_to_conditions_process.html)

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
- Sources:
    * [kratos/python_scripts/process_factory.py](https://github.com/KratosMultiphysics/Kratos/tree/master/kratos/python_scripts/process_factory.py)
    * [GitHub: ApplyConstantVectorValueProcess](https://github.com/KratosMultiphysics/Kratos/tree/master/kratos/processes/apply_constant_vectorvalue_process.cpp)
    * [Doxygen: ApplyConstantVectorValueProcess](https://kratos-docs.onrender.com/dc/d28/class_kratos_1_1_apply_constant_vector_value_process.html)