---
title: Assign Gravity
keywords: mpm gravity
tags: [mpm gravity]
sidebar: mpm_application
summary: 
---

The python process [`AssignGravityToMaterialPointProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/assign_gravity_to_material_point_process.py) assigns a vector value to the material point element variables `MP_VOLUME_ACCELERATION` and `MP_ACCELERATION` before the execution of the solution loop. As the name suggests, this process is intended to assign gravity acceleration to the material point elements of a given model part.

## Parameters & Defaults

```json
{
    "kratos_module" : "KratosMultiphysics.MPMApplication",
    "python_module" : "assign_gravity_to_material_point_process",
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
Number identifying the modulus of the imposed acceleration.

##### `direction`
Three-component vector identifying the direction of the applied acceleration.

{% include important.html content='The vector assigned to the variables `MP_ACCELERATION` and `MP_VOLUME_ACCELERATION` is given by `modulus * direction`, i.e., the scalar value `modulus` is multiplied to each component of the vector `direction`. In particular, note that the input vector `direction` is not normalized.' %}

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/assign_gravity_to_material_point_process.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/assign_gravity_to_material_point_process.py)
