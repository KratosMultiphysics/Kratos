---
title: Skin Detection
keywords: process core
tags: [skin detection process]
sidebar: kratos_core_processes
summary: 
---

# Skin Detection

## Description

This process detects the skin of a given modelpart

## Parameters & Defaults

```json
{
    "model_part_name"                       : "Main",
    "computing_model_part_name"             : "computing_domain",
    "recursive_detection"                   : false,
    "name_auxiliar_model_part"              : "SkinModelPart",
    "name_auxiliar_condition"               : "Condition",
    "list_model_parts_to_assign_conditions" : [],
    "echo_level"                            : 0
}
```

##### `model_part_name` 
Name of the source modelpart

##### `computing_model_part_name` 
Name of the auxiliar computing modelpart

##### `recursive_detection` 
if set to true will detect the skin of the modelpart and its submodelparts, if false only skin of the top level modelpart will be calculated.

##### `name_auxiliar_model_part` 
Name of the modelpart in which the skin will be calculated

##### `name_auxiliar_condition` 
Name of the condition that that will be used in the calculated skin

##### `list_model_parts_to_assign_conditions`
If provided, defines the list of modelparts in which the skin will be transfered from the `name_auxiliar_model_part`

##### `echo_level`
Output level of the process. The higher the number, the more detailed the output. (`0`, `1`, `2`)
