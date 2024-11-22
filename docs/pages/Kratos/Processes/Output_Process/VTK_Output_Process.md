---
title: 
keywords: vtk output process core
tags: [process vtk output process]
sidebar: kratos_core_processes
summary: 
---

# VTK Output process

## Description

The VTK output process specializes in generating an output in the vtk format, as its names indicates. This output is ment to be used with visualizers like paraview. 

## Execution

This process is executed in the following hooks:

##### `__init__`

Creates the folder and data structures needed for the VTK output to be generated.

##### `IsOutputStep`

Decided if the current step will produce any output

##### `PrintOutput`

Prints the output for the current step.


## Parameters & Defaults

```json
 {
    "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
    "file_format"                                 : "binary",
    "output_precision"                            : 7,
    "output_control_type"                         : "step",
    "output_interval"                             : 1.0,
    "output_sub_model_parts"                      : false,
    "output_path"                                 : "VTK_Output",
    "custom_name_prefix"                          : "",
    "custom_name_postfix"                         : "",
    "save_output_files_in_folder"                 : true,
    "entity_type"                                 : "automatic",
    "write_deformed_configuration"                : false,
    "write_ids"                                   : false,
    "nodal_solution_step_data_variables"          : [],
    "nodal_data_value_variables"                  : [],
    "nodal_flags"                                 : [],
    "element_data_value_variables"                : [],
    "element_flags"                               : [],
    "condition_data_value_variables"              : [],
    "condition_flags"                             : [],
    "gauss_point_variables_extrapolated_to_nodes" : [],
    "gauss_point_variables_in_elements"           : []
}
```

##### `model_part_name` 
Name of the modelpart in wich the process will be applied.

##### `file_format` 
Chooses the file format. Currently available options are:
- `binary`: Binary format. Smaller file sizes, impossible to read by a human 
- `ascii`: Text format. Larger file sizes, able to read it. Usefull for debug

##### `output_precision` 
Selects the number of decimals that will be printed. Maxmum number of relevant decimals is 16.

##### `output_control_type` 
Selects the control type of the output. It can be `step` or `time`

##### `output_interval` 
Selects the amout of `output_control_type`'s that need to happend before printing. (Ex 1 step, or 2.0 times)

##### `output_sub_model_parts` 
If set to `true`, will print the diefferent sub_modelparts as different outputs. If set to `false`, will print as a unified modelpart.

##### `output_path`
Path where the files will be written

##### `custom_name_prefix` 
Custom prefix that will be added to the output file: `MY_CUSTOM_PREFIX_output`

##### `custom_name_postfix` 
Custom suffix that will be added to the output file: `output_MY_CUSTOM_`

##### `save_output_files_in_folder`
If set to `true` will save all outputs of different processes inside diferent folders. Otherwise will write everything in `output_path`. usefull for MPI.

##### `entity_type` 
Selects which entity to print. Values can be:
- `Automatic` lets the output process decide
- `Element` Print Elements
- `Condition` Print Conditions

##### `write_deformed_configuration`
If set to `true` will write the deformed coordinates for the entities selected. Otherwise will print the initial coordinates.

##### `write_ids`
If set to `true` will write the id's of the entity in the `KRATOS_NODE_ID`, `KRATOS_ELEMENT_ID`, `KRATOS_CONDITION_ID` results.

##### `nodal_solution_step_data_variables` 
List of nodal variables from the historical database that will be printed.

##### `nodal_data_value_variables` 
List of nodal variables from the non-historical database that will be printed.

##### `nodal_flags` 
List of nodal flags that will be written

##### `element_data_value_variables` 
List of element variables from the non-historical database that will be printed.

##### `element_flags` 
List of element flags that will be written

##### `condition_data_value_variables` 
List of condition variables from the non-historical database that will be printed.

##### `condition_flags` 
List of condition flags that will be written

##### `gauss_point_variables_extrapolated_to_nodes` 
List of gauss point variables that will be exported as nodal results interpolating its value to the node position.

##### `gauss_point_variables_in_elements` 
List of gauss point variables that will be exported as elemental results interpolating its value to the element volume/surface.
