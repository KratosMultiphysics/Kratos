---
title: Find Global Nodal Entity Neighbours
keywords: process core
tags: [find global nodal entity neighbours process]
sidebar: kratos_core_processes
summary: 
---

# Find Global Nodal Entity Neighbours

## Description

Finds the entities that are neighbouring all the nodes from a modelpart. The result is stored in the `NEIGHBOUR_ELEMENTS` or 
`NEIGHBOUR_CONDITIONS` dependings on the entry points.

This function is exposed with the following entry points:
- `FindGlobalNodalElementalNeighboursProcess`
- `FindGlobalNodalConditionNeighboursProcess`

## Parameters & Defaults

```json
{
    "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME"
}
```

##### `model_part_name`:
Name of the modelpart in which the search will be performed