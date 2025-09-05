---
title: Calculate Distance To Path
keywords: process core
tags: [calculate distance to path process]
sidebar: kratos_core_processes
summary: 
---

# Calculate Distance To Path

## Description

This process calculates the distance from a path modelpart to a given modelpart.

## Parameters & Defaults

```json
{
    "distance_model_part_name" : "",
    "path_model_part_name"     : "",
    "distance_variable_name"   : "DISTANCE",
    "brute_force_calculation"  : false,
    "radius_path"              : 0.0,
    "distance_tolerance"       : 1.0e-9,
    "search_parameters"        :  {
        "allocation_size"         : 100,
        "bucket_size"             : 4,
        "search_factor"           : 2.0,
        "search_increment_factor" : 1.5
    }
}
```

##### `distance_model_part_name`
Name of the modelpart into wich the distance is going to be calculated.

##### `path_model_part_name`
Name of the modelpart containing the path.

##### `distance_variable_name`
Variable to store the Distance. Default `DISTANCE`

##### `brute_force_calculation`
If set to `true` will calculate the distance without using an advance search structure. Default to `false` 

##### `radius_path`
Radius of the path search. Larger values will expand the search radius.

##### `distance_tolerance`
Tolerance for the intersection between an element and a path.

##### `search_parameters`
Parameters for the Search Structure. See [SearchStructures](404.md)
