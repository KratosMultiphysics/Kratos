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
        "interval"             : [0.0, 1e30],
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

##### `interval`
Print the energy of the requested model part if the simulation time is within this interval.

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
- `"file_name"`: name of the output file.
- `"output_path"`: path of the output file.
- `"file_extension"`: extension of the output file.

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_write_energy_output_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_write_energy_output_process.py)
