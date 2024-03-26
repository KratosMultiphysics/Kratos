---
title: 
keywords: gid output process core
tags: [process gid output process]
sidebar: kratos_core_processes
summary: 
---

# GiD Output process

## Description

The GiD output process specializes in generating an output in GiD's propietary format. This output is meant to be used with GiD Post only.

## Execution

This process is executed in the following hooks:

##### `ExecuteInitialize`

Initializes the data strcutres, folders and files necessary to execute the process.

If the multifile flag is set to `false` writes the initial mesh file and leaves it open for writing.

##### `ExecuteBeforeSolutionLoop`

Writes a mesh for that particular step.
If the multifile flag is set to `ture` will write the possible placeholders for the values that may be written if its an output step

##### `IsOutputStep`

Decided if the current step will produce any output

##### `PrintOutput`

Prints the output for the current step.

##### `ExecuteFinalize`

Closes the mesh file and generate the `.lst` files containing the output names of all the files written if the flag was selected.


## Parameters & Defaults

```json
{
    "gidpost_flags": {
        "GiDPostMode": "GiD_PostBinary",
        "WriteDeformedMeshFlag": "WriteUndeformed",
        "WriteConditionsFlag": "WriteElementsOnly",
        "MultiFileFlag": "SingleFile"
    },
    "file_label": "time",
    "time_label_format": "{:.12f}",
    "output_control_type": "step",
    "output_interval": 1.0,
    "flush_after_output": false,
    "body_output": true,
    "node_output": false,
    "skin_output": false,
    "plane_output": [],
    "nodal_results": [],
    "nodal_nonhistorical_results": [],
    "nodal_flags_results": [],
    "elemental_conditional_flags_results": [],
    "gauss_point_results": [],
    "additional_list_files": []
}
```

##### `gidpost_flags` 

```json
{
    "GiDPostMode": "GiD_PostBinary",
    "WriteDeformedMeshFlag": "WriteUndeformed",
    "WriteConditionsFlag": "WriteElementsOnly",
    "MultiFileFlag": "SingleFile"
}
```

Describes the different parameters for the file output.

**`GiDPostMode`**: Chooses the file format. Currently available options are:

`GiD_PostAscii`: Text format.
`GiD_PostAsciiZipped`: Text format compressed.
`GiD_PostBinary`: Binary format.
`GiD_PostHDF5`: HDF5 compliant format.

**`WriteDeformedMeshFlag`**: Describes if the coordinates are written deformed or not:

`WriteDeformed`: Deformed
`WriteUndeformed`: Undeformed

**`WriteConditionsFlag`**: Decides if elements, conditions or both are written:

`WriteConditions`: Conditions and Elements
`WriteElementsOnly`: Only Elements
`WriteConditionsOnly`: Only Conditions

**`MultiFileFlag`**: Decides if the output will generate one single file or multiple files (one per output call)

`SingleFile`: One single file
`MultipleFiles`: One file per call

##### `file_label` 
String identifying the name of the output

##### `time_label_format` 
Prefix in the filename to indicate the output it belongs to when multiple files are written.

##### `output_control_type` 
Selects the control type of the output. It can be `step` or `time`

##### `output_interval` 
Selects the amout of `output_control_type`'s that need to happend before printing. (Ex 1 step, or 2.0 times)

##### `flush_after_output` 
If set to `true` will immediatly print the results to the file. Otherwise will let the OS manage the writting times.

##### `body_output` 
If set to true will print the values of the interal body entities.

##### `node_output` 
If set to true will print the values of the nodal entities.

##### `skin_output` 
If set to true will print the values of the interal external skin entities.

##### `plane_output` 
If set to true will print the values of the custom output plnes entities.

##### `nodal_results` 
List of nodal variables from the historical database that will be printed.

##### `nodal_nonhistorical_results` 
List of nodal variables from the non-historical database that will be printed.

##### `nodal_flags_results` 
List of nodal flags that will be written

##### `elemental_conditional_flags_results` 
List of element and condition flags that will be written

##### `gauss_point_results` 
List of gauss point variables that will be printed.

##### `additional_list_files` 
