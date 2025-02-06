---
title: Dirichlet Boundary Condition
keywords: mpm fixed displacement
tags: [mpm fixed displacement]
sidebar: mpm_application
summary: 
---

The [`ApplyMPMParticleDirichletConditionProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/apply_mpm_particle_dirichlet_condition_process.py) imposes non-conforming Dirichlet boundary conditions by means of the Penalty method, which is implemened by boundary material points conditions.

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "apply_mpm_particle_dirichlet_condition_process",
    "Parameters"    : {
        "model_part_name"                   : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "material_points_per_condition"     : 0,
        "is_equal_distributed"              : false,
        "penalty_factor"                    : 0,
        "variable_name"                     : "DISPLACEMENT",
        "constrained"                       : "fixed",
        "value"                             : [0.0, "0*t", 0.0],
        "interval"                          : [0.0, 1e30],
        "option"                            : ""
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart to which the boundary conditions will be applied.

##### `material_points_per_condition`
Defines the number of material point conditions used to replace each line or surface condition where Dirichlet boundary conditions must be imposed.

##### `is_equal_distributed`
If `true`, material point conditions are equally distributed over the original (line or surface) condition, otherwise their position is determined by Gauss quadrature nodes.

##### `penalty_factor`
Penalty coefficient used in the penalty method to enforce Dirichlet boundary conditions.

##### `variable_name`
String identifying the variable whose value is to be imposed using the Penalty method.
Admissible values are `"DISPLACEMENT"`, `"VELOCITY"` and `"ACCELERATION"`.

##### `constrained`
Specifies the type of Dirichlet boundary condition to be imposed. Admissible values are:
* `"fixed"` - the variable is fully constrained
* `"contact"` - movement is constrained only in the case of contact with the non-conforming boundary
* `"slip"` - constrains the movement along the normal direction while not restricting the others
* `"contact_slip"` - combines both contact and slip constraints

##### `value`
Value to be assigned to the variable specified by `variable_name`. Each component of the list can be either a number or a string representing a function depending on space and/or time.

##### `interval`
Time interval in which the process applies.

##### `option`
If set to `"flip_normal"`, the normal of each condition is flipped.

## Output Process File

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/apply_mpm_particle_dirichlet_condition_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/apply_mpm_particle_dirichlet_condition_process.py)
