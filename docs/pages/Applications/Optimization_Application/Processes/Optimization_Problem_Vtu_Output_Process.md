---
title: Optimization Problem Vtu Output Process
keywords: 
tags: [output, vtu, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

```OptimizationProblemVtuOutputProcess``` is used to output field data in the ```OptimizationProblem``` data container in the format vtu.

The optimization iteration data is taken from the [static data in the optimization problem](../General/Optimization_Problem.html#static-data-storage) since field data is stored in this container.

## Json settings
Following json snippet explains a single use case.
```json
{
    "type": "optimization_problem_vtu_output_process",
    "module": "KratosMultiphysics.OptimizationApplication.processes",
    "settings": {
        "file_name"                   : "<model_part_full_name>_<step>",
        "file_format"                 : "binary",
        "output_path"                 : "Optimization_Results",
        "save_output_files_in_folder" : true,
        "write_deformed_configuration": false,
        "list_of_output_components"   : ["all"],
        "output_precision"            : 7,
        "output_interval"             : 1,
        "echo_level"                  : 0
    }
}
```

| Option | Allowed values |
| ------------- | ------------- |
| type  | "optimization_problem_vtu_output_process"  |
| module  | "KratosMultiphysics.OptimizationApplication.processes" |
| file_name  | Output file name. It is allowed to use <model_part_full_name> and <step> tags. Then <model_part_full_name> tag will be replaced by the model part full name of the data field and <step> tag will be replaced by the optimization iteration.|
| file_format | "ascii" to write vtu in ascii format or "binary" to write vtu in binary format. |
| save_output_files_in_folder | Will create the output path |
| output_path | Output path to be used to store all the vtu files. |
| write_deformed_configuration | true to write the deformed configuration, false to write the initial configuration of the mesh |
| list_of_output_components | list of component names to write data. You can specify here either "all" to write data from all the components, or specify names of the response functions, controls, execution policies to write data from those components only. |
| output_precision | Precision of the data fields to be written out |
| output_interval | Interval between two consecutive vtu output writes |
| echo_level | Echo level |

## Source files
* [Doxygen](TODO) TODO
* [applications/OptimizationApplication/python_scripts/processes/optimization_problem_vtu_output_process.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/processes/optimization_problem_vtu_output_process.py)