---
title: MPM Json Output Process
keywords: mpm json output process core
tags: [mpm process json output process]
sidebar: mpm_application
summary: 
---

# MPM Json Output process

The [`MPMJsonOutputProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_json_output_process.py) exports the values of variables associated with material point **elements** in `JSON` format.

{% include important.html content='Although the `MPMJsonOutputProcess` behaves like an output process and is discussed in the output processes section, it must be included in the `"processes"` section of the `ProjectParameters.json` input file, not in the `"output_processes"` section. This requirement arises because the process extends the `KratosMultiphysics.JsonOutputProcess` class, which in turns extends `KratosMultiphysics.Process` rather than `KratosMultiphysics.OutputProcess`.' %}

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "mpm_json_output_process",
    "Parameters"    : {
        "model_part_name"               : "",
        "sub_model_part_name"           : "",
        "output_file_name"              : "",
        "gauss_points_output_variables" : [],
        "check_for_flag"                : "",
        "time_frequency"                : 1.00,
        "resultant_solution"            : false
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart.

##### `sub_model_part_name`
String identifying the name of the target SubModelPart. If not present, the whole ModelPart will be the target.

##### `output_file_name`
Name of the output file.

##### `gauss_points_output_variables`
List of variables associated to material point **elements** (e.g., `MP_MASS`, `MP_DISPLACEMENT`, `MP_VELOCITY`, ...) to be printed.

##### `check_for_flag`
A string specifying a `Flag` variable implemented in Kratos. The process
evaluates this flag for each material point element and exports the requested
data only if the flag is set to `true`. This parameter is optional; if no flag
is provided, data for all material points will be exported.

##### `time_frequency`
Defines the time interval for writing data to the file, thus determining how frequently the simulation outputs are exported.

##### `resultant_solution`
Determines how the requested variables are output for material point elements
with the `check_for_flag` condition set to `true`. If `false`, the process writes
the values of the variables separately for each material point element. If
`true`, the process calculates and writes, at each time step, the sum of the
variable values across all material points that meet the `check_for_flag`
condition.

Example. When set to `true`, the output has the form
```json
{
    "TIME" : [
        0.1,
        0.2,
        0.3
    ],
    "RESULTANT" : {
        "MP_MASS" : [
            10.0,   // sum of variable MP_MASS among all the material point elements at time 0.1
            10.0,
            10.0
        ]
    }
}
```

When set to `false`, the output has the form
```json
{
    "TIME" : [
        0.1,
        0.2,
        0.3
    ],
    "MP_1" : {
        "MP_MASS" : [
            3.0,   // value of variable MP_MASS at time 0.1 for material point elemenent with Id 1
            3.0,
            3.0
        ]
    },
    "MP_2" : {
        "MP_MASS" : [
            5.0,
            5.0,
            5.0
        ]
    },
    "MP_3" : {
        "MP_MASS" : [
            2.0,
            2.0,
            2.0
        ]
    }
}
```

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_json_output_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_json_output_process.py)
