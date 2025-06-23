---
title: Check Same Modelpart Using Skin Distance
keywords: process core
tags: [check same modelpart using skin distance process]
sidebar: kratos_core_processes
summary: 
---

# Check Same Modelpart Using Skin Distance

## Description

Checks that two modelparts are the same using the skin distance of both modelparts over an auxuliary background mesh.

Please note that checks process do not produce an output, and will only print its result.

## Parameters & Defaults

```json
{
    "skin_model_part_1_name"              : "PLEASE_SPECIFY_SKIN_MODEL_PART_2_NAME",
    "skin_model_part_2_name"              : "PLEASE_SPECIFY_SKIN_MODEL_PART_2_NAME",
    "tolerance"                           : 1.0e-3,
    "bounding_box_scale_factor"           : 1.5,
    "number_of_divisions_background_mesh" : 30,
    "discontinuous_distance_settings": {}
}
```

##### `skin_model_part_1_name` 
Name of the fist modelpart.

##### `skin_model_part_1_name` 
Name of the second modelpart.

##### `tolerance` 
Tolerance of the mesh similarity.

##### `bounding_box_scale_factor` 
Scale factor of the boundingbox containing the auxiliar background mesh.

##### `number_of_divisions_background_mesh` 
Number of divisions of the background mesh.
