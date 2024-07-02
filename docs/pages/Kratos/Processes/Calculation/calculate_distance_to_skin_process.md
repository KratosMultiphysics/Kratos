---
title: Calculate Distance To Skin
keywords: process core
tags: [calculate distance to skin process]
sidebar: kratos_core_processes
summary: 
---

# Calculate Distance To Skin

## Description

**Warning**: This process does not follow the standard Process interface.

This process takes a volume model part (with tetrahedra mesh) and a skin model part (with triangle mesh) and
and calculates the distance to the skin for all the elements and nodes of the volume model part using elemental discontinuous distances..

## Parameters & Defaults

```json
{
    "distance_variable"              : "DISTANCE",
    "distance_database"              : "nodal_historical",
    "ray_casting_relative_tolerance" : 1.0e-8
}
```

##### `distance_variable`
Variable to store the Distance. Default `DISTANCE`

##### `distance_database`
database in which the distance will be stored (`"nodal_historical"`, `"nodal_non_historical"`).  

##### `ray_casting_relative_tolerance`
Relative tolerance for the intersection of the rays using during the raycasting phase.
