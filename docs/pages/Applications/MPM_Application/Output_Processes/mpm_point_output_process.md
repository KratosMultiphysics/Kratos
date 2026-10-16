---
title: MPM Point Output Process
keywords: mpm point output process core
tags: [mpm process point output process]
sidebar: mpm_application
summary: 
---

# MPM Point Output process

The [`MPMPointOutputProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_point_output_process.py) writes results from a material point to a file. The output can be requested for material point elements and material point conditions.

## Parameters & Defaults

```json
 {
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "mpm_point_output_process",
    "Parameters"    : {
        "model_part_name"      : "",
        "entity_type"          : "element",
        "output_interval"      : 1.0,
        "output_control_type"  : "step",
        "position"             : [],
        "output_variables"     : [],
        "search_tolerance"     : 1e-6,
        "print_format"         : ".8f",
        "output_file_settings" : {}
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart.

##### `entity_type`
Selects which entity to search for. Values can be:
- `"element"`: Search material point element.
- `"condition"`: Search material point condition.

##### `output_control_type`
Determines how the output is controlled and printed during the simulation. Acceptable values are:
* `step`: the output is generated at regular intervals, specified as a fixed number of time steps;
* `time`: the output is generated at regular intervals, specified in seconds.

Choose `step` for time-step-based output control, or `time` for time-based output control.

##### `output_interval`
Selects the amount of `output_control_type` that needs to happen before printing.

Example. The following parameters are used for printing the output at each time step
```json
"output_control_type" : "step",
"output_interval"     : 1,
```

Example. The following parameters are used for printing the output every 2.0s of the simulation
```json
"output_control_type" : "time",
"output_interval"     : 2.0,
```

##### `position`
A list of 3 coordinates representing the position of the material point in space.

##### `output_variables`
List of variables to be printed to the file. Admissible material point *element* variables
are listed [here](../Material_Point_Variables/element_variables), whereas admissible material point *condition* variables are listed [here](../Material_Point_Variables/condition_variables).

##### `search_tolerance`
Tolerance value used during the search to find the material point (element or condition) closest to the specified `position`.

##### `print_format`
Format string used for printing numerical values to the output file (e.g., `".8f"` for 8 decimal places).

##### `output_file_settings`
Settings for the output file handling. This is managed by the Kratos [`TimeBasedAsciiFileWriterUtility`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/time_based_ascii_file_writer_utility.py):
```json
{
    "file_name"         : "",
    "output_path"       : "",
    "file_extension"    : "dat"
}
```
- `"file_name"`: name of the output file;
- `"output_path"`: path of the output file;
- `"file_extension"`: extension of the output file.

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_point_output_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_point_output_process.py)
