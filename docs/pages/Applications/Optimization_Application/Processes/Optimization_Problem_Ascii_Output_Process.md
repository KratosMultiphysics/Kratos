---
title: Optimization Problem Ascii Output Process
keywords: 
tags: [output, ascii, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

```OptimizationProblemAsciiOutputProcess`` is used to output ASCII data in the ```OptimizationProblem``` data container to a single file, where each optimization iteration data will be written to a new line in the output file.

The optimization iteration data is taken from the [historical data storage of the optimization problem](../General/Optimization_Problem.html#historical-data-storage). The [static data in the optimization problem](../General/Optimization_Problem.html#static-data-storage) is used as the initial data when writing.

The output will be written in the CSV format.

## Json settings
Following json snippet explains a single use case.
```json
{
    "type": "optimization_problem_ascii_output_process",
    "module": "KratosMultiphysics.OptimizationApplication.processes",
    "settings": {
        "output_file_name"         : "SPECIFY_OUTPUT_FILE_NAME",
        "write_kratos_version"     : true,
        "write_time_stamp"         : true,
        "write_initial_values"     : true,
        "list_of_output_components": ["all"],
        "format_info": {
            "int_length"     : 7,
            "float_precision": 9,
            "bool_values"    : ["no", "yes"],
            "string_length"  : 10
        }
    }
}
```

| Option | Allowed values |
| ------------- | ------------- |
| type  | "optimization_problem_ascii_output_process"  |
| module  | "KratosMultiphysics.OptimizationApplication.processes" |
| output_file_name  | Output file name. |
| write_kratos_version | true to write the Kratos version used in the optimization analysis, false will not write anything. |
| write_time_stamp | true to write the time stamp of the starting point, false will not write anything |
| write_initial_values | true to write the initial data. |
| list_of_output_components | list of component names to write data. You can specify here either "all" to write data from all the components, or specify names of the response functions, controls, execution policies to write data from those components only. |
| int_length | Length of the integer values when converted to strings |
| float_precision | Floating point value precision to be used when converting to strings |
| bool_values | The keywords to be used if bool is true and false. |
| string_length | String length to be used to write string data. |

## Source files
* [Doxygen](TODO) TODO
* [applications/OptimizationApplication/python_scripts/processes/optimization_problem_ascii_output_process.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/processes/optimization_problem_ascii_output_process.py)


