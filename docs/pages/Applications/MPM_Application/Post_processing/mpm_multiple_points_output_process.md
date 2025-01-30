---
title: MPM Multiple Points Output Process
keywords: mpm point points output process core
tags: [mpm point points output process]
sidebar: mpm_application
summary: 
---

# MPM Multiple Points Output Process

The [`MPMMultiplePointsOutputProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_multiple_points_output_process.py) writes variables associated to **multiple** material points (element or condition) in the model to a file. Each material point (element or condition) is identified by searching the entity closest to the specified location, provided it lies within a defined tolerance. It is based on the [`MPMPointOutputProcess`](./mpm_point_output_process).

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "mpm_multiple_points_output_process",
    "Parameters"    : {
        "model_part_name"       : "",
        "entity_type"           : "element|condition",
        "interval"              : [0.0, 1e30],
        "positions"             : [[]],
        "search_tolerance"      : 1e-6,
        "output_variables"      : [],
        "print_format"          : "",
        "output_file_settings"  : {
            "file_name"         : "",
            "output_path"       : "",
            "write_buffer_size" : -1,
            "file_extension"    : "dat"
        }
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart.

##### `entity_type`
String identifying the entity to search. Available options are `element` and `condition`.

##### `interval`
Print the requested variables if the simulation time is within this interval.

##### `output_variables`
Selects the variables to print. For every scalar variable (e.g., `MP_DENSITY`) a column will be printed, whereas for every array variable (e.g., `MP_DISPLACEMENT`) a column for every component will be printed (i.e., `MP_DISPLACEMENT_X`, `MP_DISPLACEMENT_Y` and `MP_DISPLACEMENT_Z`).

##### `positions`
Array of the positions of the material point entities (elements or conditions) with respect to the initial configuration.

##### `search_tolerance`
Tolerance of the search. For each entry of the array `positions`, we select the entity that is inside a sphere of radius given by `tolerance` and centered on the given position.
If several entities meet this condition, we select the one that is closest to the given position.

##### `print_format`
Overrides the default print format. Must be a string with a token placeholder `{}`.

##### `print_output_file_settings`
Defines the settings for the `TimeBasesAsciiFileWriterUtility`.

```json
{
    "file_name"         : "",
    "output_path"       : "",
    "file_extension"    : "dat"
}
```

- `file_name`: name of the output file.
- `output_path`: path of the output file.
- `file_extension`: extension of the output file.

## Output Process File

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_multiple_points_output_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_multiple_points_output_process.py)
