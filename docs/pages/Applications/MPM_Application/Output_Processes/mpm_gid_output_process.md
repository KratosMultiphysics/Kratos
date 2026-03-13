---
title: MPM GiD Output Process
keywords: mpm gid output process core
tags: [mpm process gid output process]
sidebar: mpm_application
summary: 
---

# MPM GiD Output process

The [`MPMGiDOutputProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_gid_output_process.py) generates an output in `GiD`'s propietary format.

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "mpm_gid_output_process",
    "Parameters"    : {
        "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "postprocess_parameters" : {
            "result_file_configuration" : {
                "output_control_type" : "step",
                "output_interval"     : 1,
                "gauss_point_results" : []
            }
        },
        "output_name" : ""
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

##### `gauss_point_results`
List of variables associated to material points (e.g., `MP_DISPLACEMENT`, `MP_VELOCITY`, ...) to be printed.

##### `output_name`
Name of the output file.

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_gid_output_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_gid_output_process.py)
