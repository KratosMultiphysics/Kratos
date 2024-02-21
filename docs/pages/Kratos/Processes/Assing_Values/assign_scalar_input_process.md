---
title: Assign Scalar Input
keywords: process core
tags: [assign scalar input to entities process]
sidebar: kratos_core_processes
summary: 
---

# Assign Scalar Input

## Description

This process assigns a value from an input to a variable belonging to all of the entities in a given mesh

## Execution

This process is executed in the following hooks:

#### `ExecuteInitializeSolutionStep`

Assigns the value.

## Parameters & Defaults

```json
{
    "model_part_name"    : "please_specify_model_part_name",
    "mesh_id"            : 0,
    "variable_name"      : "SPECIFY_VARIABLE_NAME",
    "interval"           : ["begin", "end"],
    "file"               : "",
    "historical"         : false,
    "transfer_algorithm" : "nearest_neighbour",
    "entities"           : []
}
```

##### `model_part_name`
Name of the modelpart in wich the field variable will be applied

##### `mesh_id`
Id of the internal mesh to which the process will be applied. Default `0`.

##### `variable_name`
Name of the variable in which the field value will be applied.

##### `interval`
Interval of time in which the process will be applied.

##### `file`
Input file with the input values. currently accepts `json`, `txt` and `csv`

##### `historical`
Selects if the value is applied to the historical (`true`) or the non-historical (`false`) databases. Default `false`

##### `transfer_algorithm`
Searh algortihm. Curently only accepts `nearest_neighbour`.

##### `entities`
List of entities into which the value will be applies. Accepts: `nodes`, `elements`, `conditions`

## File Formats

### JSON

A JSON file containing 
- A field `TIME` with an array of time steps in which the value will be applied.
- A list of fields containing:
    - `COORDINATES` with the coordinates into which the value will be applied
    - `VALUES` a field containing. Please notice that while applying from json, is possible to define multiple values which will override the `variable_name` definided in the properties.
        - The name of the variable to which will be applied.
        - A list of values. Needs to be of the same size as the `TIME` field in the root. 

#### Example

```json
{
    "TIME" : [0.1,0.2,0.3,0.4,0.5],
    "1"    : {
        "COORDINATES" : [0.75, 0.75, 0.0],
        "VALUES"      : { "PRESSURE" : [1.0,2.0,4.0,8.0,16.0] }
    },
    "2"    : {
        "COORDINATES" : [0.25, 0.25, 0.0],
        "VALUES"      : { "PRESSURE" : [2.0,4.0,8.0,16.0,32.0] }
    }
}

```

### TXT

A text file containing:
- A header lines with the `time` keyword and a series of coordinates that identify the position in which the value will be applied.
- List of rows containing the time step of application and the values. The list of values needs to be of the same size as the number of coordinates in the header.

#### Example

```txt
time 	(0.75, 0.75, 0.0) 	(0.25, 0.25, 0.0)
0.1 	1.0 	2.0
0.2 	2.0 	4.0
0.3 	4.0 	8.0
0.4 	8.0 	16.0
0.5  	16.0 	32.0
```

### CSV

Header containing `#time velocity` and a list of values

#### Example

```csv
#time velocity
0,  0
1,  0.1
2,  0.3
```
