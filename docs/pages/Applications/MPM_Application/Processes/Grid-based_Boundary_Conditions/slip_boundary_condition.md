---
title: Slip
keywords: mpm slip
tags: [mpm slip]
sidebar: mpm_application
summary: 
---

The [`ApplyMPMSlipBoundaryProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/apply_mpm_slip_boundary_process.py) imposes zero displacement in the direction orthogonal to the conditions of a background grid submodel part (no constraints in the tangent direction).

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "apply_mpm_slip_boundary_process",
    "Parameters"    : {
        "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "option"                    : "",
        "avoid_recomputing_normals" : true
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart.

##### `option`
If set to `"flip_normal"`, the normal of each node of the target model part is flipped.

##### `avoid_recomputing_normals`
If `false`, the process recomputes the normals of each node of the target model part in the initialization of every solution step.

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/apply_mpm_slip_boundary_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/apply_mpm_slip_boundary_process.py)
