---
title: Point Load Condition
keywords: mpm point load
tags: [mpm point load]
sidebar: mpm_application
summary: 
---

The [`ApplyMPMParticleNeumannConditionProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/apply_mpm_particle_neumann_condition_process.py) imposes a point load on a (moving) material point condition.

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "apply_mpm_particle_neumann_condition_process",
    "Parameters"    : {
        "model_part_name"   : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "variable_name"     : "POINT_LOAD",
        "constrained"       : "fixed",
        "modulus"           : 1.0,
        "direction"         : [0.0, 0.0, 0.0],
        "interval"          : [0.0, 1e30],
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart to which the process is applied.

##### `variable_name`
String identifying the variable whose value is to be imposed using the Penalty method.
The only admissible values is `POINT_LOAD`.

##### `constrained`
Admissible values are:
* `"fixed"`
* `"normal"`

##### `modulus`
Modulus of the vector to be assigned to the variable specified by `variable_name`.
The value can be either a number or a string representing a function depending on
space and/or time.

##### `direction`
Direction of the vector to be assigned to the variable specified by `variable_name`.
Each component of the list can be either a number or a string representing a function depending on space and/or time.

##### `interval`
Time interval in which the process applies.

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/apply_mpm_particle_neumann_condition_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/apply_mpm_particle_neumann_condition_process.py)
