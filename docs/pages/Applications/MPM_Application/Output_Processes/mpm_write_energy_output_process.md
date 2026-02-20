---
title: MPM Write Energy Output Process
keywords: mpm vtk output process core
tags: [mpm process vtk output process]
sidebar: mpm_application
summary: 
---

# MPM Write Energy Output process

The [`MPMWriteEnergyOutputProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_write_energy_output_process.py) generates an output file in which the **kinetic**, **potential**, **strain** and **total energy** of a given model part are written at each time step.

## Parameters & Defaults

```json
 {
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "mpm_write_energy_output_process",
    "Parameters"    : {
        "model_part_name"      : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "output_control_type"  : "step",
        "output_interval"      : 1.0,
        "print_format"         : ".8f",
        "output_file_settings" : {
            "file_name"         : "",
            "output_path"       : "",
            "file_extension"    : "dat"
        }
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart.

##### `output_control_type`
Determines how the output is controlled and printed during the simulation. Acceptable values are:
* `step`: the output is generated at regular intervals, specified as a fixed number of time steps;
* `time`: the output is generated at regular intervals, specified in seconds.

Choose `step` for time-step-based output control, or `time` for time-based output control.

##### `output_interval`
Selects the amount of `output_control_type` that needs to happend before printing.

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
String identifying the format to be used for printing the energy.

##### `output_file_settings`
Defines the settings for the [`TimeBasesAsciiFileWriterUtility`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/time_based_ascii_file_writer_utility.py):
```json
{
    "file_name"         : "",
    "output_path"       : "",
    "file_extension"    : "dat"
}
```
- `"file_name"`: name of the output file. If the name is not specified, the file name is `<model_part_name>_energy.dat`.
- `"output_path"`: path of the output file.
- `"file_extension"`: extension of the output file.

## Example
Ouput example:

```
# Energy model part 'MPMModelPart'
# time potential_energy kinetic_energy strain_energy total_energy
0.0 0.00000000 0.14000000 4.55000000 4.69000000
0.2 0.41700000 1.27500000 13.21320000 14.90520000
0.4 1.20195000 7.48890000 25.98960000 34.68045000
0.6 2.70489915 27.29754000 43.53440000 73.53683915
```

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_write_energy_output_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_write_energy_output_process.py)
