---
title: 
keywords: Point output process core
tags: [point output process]
sidebar: kratos_core_processes
summary: 
---

# Point Output Process

## Description

This process writes results from a geometrical position (point) in the model to a file. It first searches the entity containing the requested output location and then interpolates the requested variable(s). The output can be requested for elements, conditions and nodes. For nodes no geometrical interpolation is performed, the exact coordinates have to be specified. This process works in MPI as well as with restarts. It can serve as a basis for other processes (e.g. MultiplePointsOutputProcess). Furthermore it can be used for testing in MPI where the node numbers can change

## Execution

This process is executed in the following hooks:

##### `ExecuteInitialize`

Initializes the search structures, the internal search strategy (Normal, Restarted, MPI).

##### `ExecuteInitializeSolutionStep`

Search for the entities to be printed. This operation is done for every time step loop in case the entities have changed position.

##### `IsOutputStep`

Decides if the current step will produce any output

##### `PrintOutput`

Prints the output for the current step.

##### `ExecuteFinalize`

Closes the output file and finalizes the process

## Parameters & Defaults

```json
{
    "model_part_name"      : "",
    "entity_type"          : "[node|element|condition]",
    "interval"             : [0.0, 1e30],
    "position"             : [],
    "output_variables"     : [],
    "historical_value"     : true,
    "search_configuration" : "[initial|current]",
    "search_tolerance"     : 1e-6,
    "print_format"         : "",
    "output_file_settings" : {
        "file_name"  : "",
        "output_path": "",
        "write_buffer_size" : -1,
        "file_extension" : "dat"
    }
}
```


##### `model_part_name` 
String identifying the name of the target modelpart.

##### `entity_type` 
Selects the type of entity that will print the data. If the entity its a node, the process will search for the closest node from the given coordinates in `position`. Otherwise the process will perform an interpolation of the selected values in the selected coordinates over all the nodes belonging  to the entity in the selected coordinates.

##### `interval` 
Selects the control type of the output. It can be `step` or `time`

##### `position` 
The position to be printed.

##### `output_variables` 
Selects the variables to print. For every scalar variable a column will be printed, for every array variable a column for every component will be printed. This behaviour can be changed passing a `print_format` string.

##### `historical_value` 
Selects if the values printed will be taken from the historical (`true`) or the non-historical (`false`) databases. Note that selecting non-historical values form non-node elements, will still interpolate their non-historical nodal values.

##### `search_configuration` 
Selects if the search will be performed over the initial position of the entities (`initial`) or over the current position of the entities (`current`)

##### `search_tolerance` 
Tolerance of the search.

##### `print_format` 
Overrides the default print format of the entity values selected. The process will honour any string with a token placeholder `{}` but will fail if the amount of tokens in not the expected (1 for every scalar variable, 3 for every vector variable)

##### `print_output_file_settingsformat` 
See also [TimeBasedAsciiFileWriterUtility](../Utilities/time_based_ascii_file_writer_utility.md)

Defines the settings for the `TimeBasesAsciiFileWriterUtility`.

```json
{
    "file_name"  : "",
    "output_path": "",
    "write_buffer_size" : -1,
    "file_extension" : "dat"
}
```

- **`file_name`**: Name of the output file.

- **`output_path`**: Path of the output file.

- **`write_buffer_size`**: Size of the internal buffer. `-1` for system managed, and any value `>0` for a custom buffer size.

- **`file_extension`**: Extension of the output file.
