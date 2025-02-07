---
title: MPM Point Output Process
keywords: mpm Point output process core
tags: [mpm point output process]
sidebar: mpm_application
summary: 
---

# MPM Point Output Process

The [`MPMPointOutputProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_point_output_process.py) writes variables associated with a specific material point (element or condition) in the model to a file. The material point is identified by searching the entity closest to the specified location, provided it lies within a defined tolerance.

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "mpm_point_output_process",
    "Parameters"    : {
        "model_part_name"       : "",
        "entity_type"           : "element|condition",
        "interval"              : [0.0, 1e30],
        "output_variables"      : [],
        "position"              : [],
        "search_tolerance"      : 1e-6,
        "print_format"          : "",
        "output_file_settings"  : {
            "file_name"         : "",
            "output_path"       : "",
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

##### `position`
The position of the material point entity (element or condition) with respect to the initial configuration.

##### `search_tolerance`
Tolerance of the search. We select the entity that is inside a sphere of radius given by `tolerance` and centered on `position`.
If several entities meet this condition, we select the one that is closest to the given `position`.

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

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_point_output_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_point_output_process.py)
