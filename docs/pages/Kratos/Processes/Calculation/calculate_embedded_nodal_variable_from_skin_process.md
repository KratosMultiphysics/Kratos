---
title: Calculate Embedded Nodal Variable From Skin
keywords: process core
tags: [calculate embedded nodal variable from skin process]
sidebar: kratos_core_processes
summary: 
---

# Calculate Embedded Nodal Variable From Skin

## Description

This process tramsfers the value from a source variable of an skin modelpart and embbeds it into a target variable of a target modelpart.

This process needs an auxliar modelpart to make internal calculations.

## Parameters & Defaults

```json
{
    "echo_level" : 0,
    "base_model_part_name": "",
    "skin_model_part_name": "",
    "skin_variable_name": "",
    "embedded_nodal_variable_name": "",
    "buffer_position": 0,
    "aux_model_part_name": "IntersectedElementsModelPart",
    "gradient_penalty_coefficient": 0.0,
    "linear_solver_settings": {
        "preconditioner_type": "amg",
        "solver_type": "amgcl",
        "smoother_type": "ilu0",
        "krylov_type": "cg",
        "max_iteration": 1000,
        "verbosity": 0,
        "tolerance": 1e-8,
        "scaling": false,
        "block_size": 1,
        "use_block_matrices_if_possible": true
    }
}
```

##### `echo_level`
Output level of the process. The higher the number, the more detailed the output. (`0`, `1`, `2`)

##### `base_model_part_name`
Target modelpart into wich the variable value will be transfered.

##### `skin_model_part_name`
Source modelpart form which the variable will be sourced.

##### `skin_variable_name`
Name of the skin variable to be transfer.

##### `embedded_nodal_variable_name`
Name of the target variable from the target modelpart.

##### `buffer_position`
Time buffer from which the variables will be taken. Defaults to `0` (current time step)

##### `aux_model_part_name`
Name of the auxiliar modelpart

##### `gradient_penalty_coefficient`
Applies a gradient penalty coeficient to the auxiliar modelpart

##### `linear_solver_settings`
Settings for the linear solver. See [Linear Solvers]() for more info.


