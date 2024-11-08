---
title: Calculate Discontinuous Distance To Skin
keywords: process core
tags: [calculate discontinuous distance to skin process]
sidebar: kratos_core_processes
summary: 
---

# Calculate Discontinuous Distance To Skin

## Description

**Warning**: This process does not follow the standard Process interface.

This process takes a volume model part (with tetrahedra mesh) and a skin model part (with triangle mesh) and
and calcualtes the distance to the skin for all the elements and nodes of the volume model part.

## Parameters & Defaults

```json
{
    "elemental_distances_variable"                   : "ELEMENTAL_DISTANCES",
    "elemental_edge_distances_variable"              : "ELEMENTAL_EDGE_DISTANCES",
    "elemental_edge_distances_extrapolated_variable" : "ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED",
    "embedded_velocity_variable"                     : "EMBEDDED_VELOCITY",
    "calculate_elemental_edge_distances"             : false,
    "calculate_elemental_edge_distances_extrapolated": false,
    "use_positive_epsilon_for_zero_values"           : true
}
```

##### `elemental_distances_variable` 
Variable that will be use to store the distance between elements.

##### `elemental_edge_distances_variable` 
Variable that will be use to store the cut edge ratios.

##### `elemental_edge_distances_extrapolated_variable` 
Variable that will be use to store the cut edge ratios of extrapolated geometry.

##### `embedded_velocity_variable`
Variable that holds the embedded velocity of the mesh.

##### `calculate_elemental_edge_distances`
If set to `true`, will calculate the cut edge ratios for intersected elements and store them in `elemental_edge_distances_variable` 

##### `calculate_elemental_edge_distances_extrapolated`
If set to `true`, will calculate the cut edge ratios of the extrapolated geometry for intersected elements and store them in `elemental_edge_distances_extrapolated_variable` 

##### `use_positive_epsilon_for_zero_values`
If set to `true` will use positive `std::numeric_limits<double>::epsilon()` for the tolerace of the distance calculation, otherwise will use `-std::numeric_limits<double>::epsilon()`