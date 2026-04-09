---
title: MPM Grid Conforming Reaction Output Process
keywords: mpm grid-conforming reaction output process
tags: [mpm process reaction output process]
sidebar: mpm_application
summary: 
---

# MPM Grid Conforming Reaction Output process

The [`MPMGridConformingReactionOutputProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_grid_conforming_reaction_output_process.py) generates a text file containing the total reaction forces over time on **grid-conforming boundaries**.

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "mpm_grid_conforming_reaction_output_process",
    "Parameters"    : {
        "model_part_name"           : "",
        "output_interval"           : 1.0,
        "output_control_type"       : "step",
        "print_format"              : ".8f",
        "output_file_settings"      : {}
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart. **Note**: For this specific process, the model part name must start with `Background_Grid`.

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

##### `print_format`
Formatting string to specify the precision of the output values (e.g., `".8f"` for floating-point with 8 decimal places).

##### `output_file_settings`
Defines the settings for the [`TimeBasesAsciiFileWriterUtility`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/time_based_ascii_file_writer_utility.py):
```json
{
    "file_name"         : "",
    "output_path"       : "",
    "file_extension"    : "dat"
}
```
- `"file_name"`: name of the output file. If the name is not specified, the file name is `<model_part_name>_reaction`.
- `"output_path"`: path of the output file.
- `"file_extension"`: extension of the output file.

## Example
Output example:

```
# Grid conforming reaction for model part Background_Grid.grid_conforming_reaction
# Time Fx Fy Fz
0.0 0.00000000 0.00000000 0.00000000
0.2 3.80000000 7.60000000 15.20000000
0.4 5.60000000 11.20000000 22.40000000
0.6 7.40000000 14.80000000 29.60000000
```

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_grid_conforming_reaction_output_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_grid_conforming_reaction_output_process.py)
