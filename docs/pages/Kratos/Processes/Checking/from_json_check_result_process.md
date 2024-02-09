---
title: From Json Check Result
keywords: process core
tags: [from json check result process]
sidebar: kratos_core_processes
summary: 
---

# From Json Check Result

## Description

This process checks results using a json file containing the solution af given model part with a certain frequency

## Parameters & Defaults

```json
{
    "check_variables"              : [],
    "gauss_points_check_variables" : [],
    "input_file_name"              : "",
    "model_part_name"              : "",
    "sub_model_part_name"          : "",
    "check_for_flag"               : "",
    "historical_value"             : true,
    "tolerance"                    : 1e-3,
    "relative_tolerance"           : 1e-6,
    "time_frequency"               : 1.00,
    "use_node_coordinates"         : false,
    "check_only_local_entities"    : false
}
```

##### `check_variables` 
List of variables to be check.

##### `gauss_points_check_variables` 
List of gauss points to be check.

##### `input_file_name` 
Name of the file with the reference values.

##### `model_part_name` 
Name of the modelpart to be check.

##### `sub_model_part_name` 
Name of the submodelpart to be check (optional).

##### `check_for_flag` 
Name of the flag to be check. Please note that due to restrictions in the falg implementation only one flag can be checked at a time.

##### `historical_value` 
If `true` historical variables will be checked, if `false` non historical variables will be checked.

##### `tolerance` 
Absolute tolerance for the value comparison.

##### `relative_tolerance` 
Rleative tolerance for the value comparison.

##### `time_frequency` 
Amount of time between checks. If the process is called before time_frequency has passed the check will no be performed.

##### `use_node_coordinates` 
Currently unused.

##### `check_only_local_entities` 
If `true` will check only entities in the main process, otherwise will check entitites in the MPI interfaces.