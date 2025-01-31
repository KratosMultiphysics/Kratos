---
title: Assign Initial Velocity
keywords: mpm velocity
tags: [mpm velocity]
sidebar: mpm_application
summary: 
---

The python process [`AssignInitialVelocityToMaterialPointProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/assign_initial_velocity_to_material_point_process.py) assigns a vector value to the variable `MP_VELOCITY` of all the elements of a given model part before the execution of the solution loop.

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "assign_initial_velocity_to_material_point_process",
    "Parameters"    : {
        "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "modulus"         : 1.0,
        "direction"       : [0.0, 0.0, 0.0]
    }
}
```

##### `model_part_name`
String identifying the name of the target ModelPart.

##### `modulus`
Number identifying the modulus of the imposed velocity.

##### `direction`
Three-component vector identifying the direction of the applied velocity.

{% include important.html content='The vector assigned to variable the `MP_VELOCITY` is given by `modulus * direction`, i.e., the scalar value `modulus` is multiplied to each component of the vector `direction`. In particular, note that the input vector `direction` is not normalized.' %}

## Process File

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/assign_initial_velocity_to_material_point_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/assign_initial_velocity_to_material_point_process.py)
