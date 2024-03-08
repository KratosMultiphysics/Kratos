---
title: Set-up Example
keywords: json
tags: [JSON_Structure.md]
sidebar: structural_mechanics_application
summary:
---
# Overview
This page will guide the users how to configure the json for the Structural Mechanics Application.

## Problem Data
The json file starts with the problem data. It consists of the following parameters.

```json
{
    "problem_data":
    {
        "problem_name"    : "<problem_name>",
        "parallel_type"   : "OpenMP",
        "start_time"      : 0.0,
        "end_time"        : 20.0,
        "echo_level"      : 0
    },
```

## Solver Settings

### solver_type

Users can choose the solver type based on the art of the problem. The available names for the solvers are:
- static
- dynamic
- eigen_value

### domain_size

The domain size is the number of dimensions for which the problem is designed eg. for 3D problems, the domain size should be set to 3.

### model_import_settings

```json
"model_import_settings"              : {
            "input_type"       : "mdpa",
            "input_filename"   : "<path_to_your_file/file_name>"
        },
```
Import your model using the above flag. `input_type` sets the file type eg. mdpa and the `input_filename` sets the path to the input file.

### material_import_settings

```json
"material_import_settings"           : {
            "materials_filename" : "<path_to_your_file/file_name>.json"
        },
```
Import the materials parameter `json` file using the above flag.

### analysis_type

Use this flag to denote the linearity of the problem, i.e., the available options are:
- linear
- non-linear

(Not applicable for `eigen_value` solver type.)

### Example

```json
"solver_settings"          : {
        "solver_type"                        : "Static",
        "echo_level"                         : 0,
        "model_part_name" : "Structure",
        "domain_size"     : 3,
        "time_stepping"                      : {
            "time_step" : 0.4
        },
        "analysis_type"                      : "non_linear",
        "model_import_settings"              : {
            "input_type"     : "mdpa",
            "input_filename" : "truss_test/linear_3D2NTruss_plastic_compression_test"
        },
        "material_import_settings" :{
            "materials_filename": "truss_test/linear_3D2NTruss_plastic_compression_test_materials.json"
        },
        "line_search"                        : false,
        "convergence_criterion"              : "residual_criterion",
        "displacement_relative_tolerance"    : 0.0001,
        "displacement_absolute_tolerance"    : 1e-9,
        "residual_relative_tolerance"        : 1e-9,
        "residual_absolute_tolerance"        : 1e-9,
        "max_iteration"                      : 40,
        "rotation_dofs"                      : false
    },
```

## Processes

This section in the JSON configuration defines a set of operations or tasks to be performed during a simulation. These processes collectively control the behavior of the simulation, imposing constraints, applying loads, and performing result checks.

### constraints_process_list

The flags inside this list specify details about the constraints imposed on vector variables in a structural simulation.

##### Example
```json
"processes" : {
    "constraints_process_list" : [{
        "python_module" : "assign_vector_variable_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"          : "AssignVectorVariableProcess",
        "Parameters"            : {
            "mesh_id"         : 0,
            "model_part_name" : "Structure.DISPLACEMENT_dirichletXYZ",
            "variable_name"   : "DISPLACEMENT",
            "constrained"	  : [true,true,true],
            "value"           : [0.0,0.0,0.0]
        }
    },
```
Here, `Parameters` flag is a dictionary containing settings for each process in the "constraints_process_list," specifying details such as the mesh, model part, variable, constrained components, fixed values, and the time interval for constraint application.

### loads_process_list

This list configures processes that apply directed loads to specific conditions in the simulation model.

##### Example

```json
"loads_process_list"       : [{
            "python_module" : "assign_vector_by_direction_to_condition_process",
            "kratos_module" : "KratosMultiphysics",
            "Parameters"    : {
                "model_part_name" : "Structure.PointLoad3D_forceX",
                "variable_name"   : "POINT_LOAD",
                "modulus"         : 0.0,
                "direction"       : [1,0.0,0.0],
                "interval"        : [0.0,"End"]
            }
        },
```
Here, `Parameters` includes details such as the model part, variable, load magnitude, direction, and the time interval during which the process is active.

### list_other_processes

All other miscellaneous processes could be listed inside this flag.

##### Example

```json
"list_other_processes": [
    {
        "python_module"   : "from_json_check_result_process",
        "kratos_module" : "KratosMultiphysics",
        "help"                  : "",
        "process_name"          : "FromJsonCheckResultProcess",
        "Parameters"            : {
            "check_variables"  : ["DISPLACEMENT_X","DISPLACEMENT_Y","REACTION_Y","REACTION_X"],
            "input_file_name"  : "beam_test/linear_3D2NBeamCr_test_results.json",
            "model_part_name"  : "Structure",
            "time_frequency"   : 0.9
        }
    }
```
Here, all the flags are same as constraint and load process lists, however the `Parameters` flag should be adapted to have different parameters inside depending upon the corresponding process.