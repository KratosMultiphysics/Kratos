---
title: Apply Periodic Boundary Condition
keywords: process core
tags: [apply periodic boundary condition process]
sidebar: kratos_core_processes
summary: 
---

# Apply Periodic Boundary Condition

## Description

**Warning**: This is a legacy process and does not follow the standard Process interface.

This process will apply a boundary condition to the modelpart is applied to.

This process is executed on the follwing hooks:
- `ExecuteInitialize`

## Parameters & Defaults

```json
{
    "variable_names":[],
    "transformation_settings":{
        "rotation_settings":{
            "center":[0,0,0],
            "axis_of_rotation":[0.0,0.0,0.0],
            "angle_degree":0.0
        },
        "translation_settings":{
            "dir_of_translation":[0.0,0.0,0.0],
            "magnitude":0.0
        }
    },
    "search_settings":{
        "max_results":100000,
        "tolerance": 1E-6
    }
}
```