---
title: MPM Vtk Output Process
keywords: mpm vtk output process core
tags: [mpm process vtk output process]
sidebar: mpm_application
summary: 
---

# MPM Vtk Output process

The [`MPMVtkOutputProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_vtk_output_process.py) generates an output in the vtk format. This output is meant to be used with softwares such as [Paraview](https://www.paraview.org/) or [VisIt](https://visit-dav.github.io/visit-website/).

## Parameters & Defaults

```json
 {
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "mpm_vtk_output_process",
    "Parameters"    : {
        "model_part_name"                   : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "file_format"                       : "binary",
        "output_precision"                  : 7,
        "output_control_type"               : "step",
        "output_interval"                   : 1.0,
        "output_sub_model_parts"            : false,
        "output_path"                       : "MPM_VTK_Output",
        "custom_name_prefix"                : "",
        "custom_name_postfix"               : "",
        "save_output_files_in_folder"       : true,
        "entity_type"                       : "automatic",
        "write_ids"                         : false,
        "element_flags"                     : [],
        "condition_flags"                   : [],
        "gauss_point_variables_in_elements" : []
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart.

##### `file_format`
Chooses the file format. Currently available options are:
- `binary` (binary format): smaller file sizes, impossible to read by a human;
- `ascii` (text format): larger file sizes, able to read it; usefull for debug.

##### `output_precision`
Selects the number of decimals that will be printed. Maxmum number of relevant decimals is 16.

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

##### `output_sub_model_parts`
If set to `true`, the process will print the diefferent SubModelParts as different output files.
If set to `false`, the process will print everything in one file.

##### `output_path`
Path where the files will be written.

##### `custom_name_prefix`
Custom prefix that will be added to the output file: `MY_CUSTOM_PREFIX_output`.

##### `custom_name_postfix`
Custom suffix that will be added to the output file: `output_MY_CUSTOM_SUFFIX`.

##### `save_output_files_in_folder`
If set to `true` will save all outputs of different MPI processes inside diferent folders. Otherwise will write everything in `output_path`.

##### `entity_type`
Selects which entity to print. Values can be:
- `Automatic`: the output process decides whether to print elements or conditions
- `Element`: print elements
- `Condition`: print conditions

##### `write_ids`
If set to `true`, the process will write the id's of the entities.

##### `element_flags`
List of element flags that will be written.

##### `condition_flags`
List of condition flags that will be written.

##### `gauss_point_variables_in_elements`
List of material point variables to be printed.

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_vtk_output_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_vtk_output_process.py)
