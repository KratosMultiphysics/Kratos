---
title: MPM Multiple Points Output Process
keywords: mpm multiple points output process core
tags: [mpm process multiple points output process]
sidebar: mpm_application
summary: 
---

# MPM Multiple Points Output process

The [`MPMMultiplePointsOutputProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_multiple_points_output_process.py) writes results from several material points to different files. Internally, it holds objects of type [`MPMPointOutputProcess`](./mpm_point_output_process).

## Parameters & Defaults

```json
 {
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "mpm_multiple_points_output_process",
    "Parameters"    : {
        "model_part_name"       : "",
        "entity_type"           : "element",
        "output_control_type"   : "step",
        "output_interval"       : 1.0,
        "positions"             : [[]],
        "output_variables"      : [],
        "search_tolerance"      : 1e-6,
        "print_format"          : ".8f",
        "output_file_settings"  : {}
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart.

##### `entity_type`
Selects which entities to search for. Values can be:
- `"element"`: Search material point elements.
- `"condition"`: Search material point conditions.

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

##### `positions`
A matrix (list of lists) containing the 3 coordinates representing the position of each material point in space. For example, to output two points: `[[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]`.

##### `output_variables`
List of variables to be printed to the file. Admissible material point *element* variables
are listed [here](../Material_Point_Variables/element_variables), whereas admissible material point *condition* variables are listed [here](../Material_Point_Variables/condition_variables).

##### `search_tolerance`
Tolerance value used during the search to find the material points (elements or conditions) closest to the specified `positions`.

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
- `"file_name"`: base name of the output files. An index will be appended for each point (e.g., `_1`, `_2`);
- `"output_path"`: path of the output files;
- `"file_extension"`: extension of the output files.

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_multiple_points_output_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_multiple_points_output_process.py)
