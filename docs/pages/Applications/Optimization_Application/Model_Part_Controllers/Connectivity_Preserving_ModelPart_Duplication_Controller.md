---
title: Connectivity Preserving ModelPart Duplication Controller
keywords: 
tags: [connectivity preserving modeller, model part controller, optimization]
sidebar: optimization_application
summary: 
---
## Introduction

This ```ConnectivityPreservingModelPartDuplicationController``` is used to create a new model part or fill an existing model part with specified elements and conditions. All of these elements and conditions in the destination model part will have the same shared geometry and nodes with the source model part. Hence, the connectivities are preserved in the destination model part.

## Json settings
Following json-snippet illustrates an example use case
```json
{
    "type": "connectivity_preserving_model_part_duplication_controller",
    "module": "KratosMultiphysics.OptimizationApplication.model_part_controllers",
    "settings": {
        "source_model_part_name"      : "AdjointStructure",
        "destination_model_part_name" : "Structure",
        "destination_element_name"    : "SolidElement3D4N",
        "destination_condition_name"  : "SurfaceCondition3D3N"
    }
}
```

| Option | Allowed values |
| ------------- | ------------- |
| type  | "connectivity_preserving_model_part_duplication_controller"  |
| module  | "KratosMultiphysics.OptimizationApplication.model_part_controllers"  |
| source_model_part_name  | Source model part name to read the nodes and the connectivities. |
| destination_model_part_name  | Destination model part name. If not found, it will be created.|
| destination_element_name  | Elements to be used in the destination model part |
| destination_condition_name  | Conditions to be used in the destination model part |

## Source files
* [Doxygen](TODO) TODO
* [applications/OptimizationApplication/python_scripts/model_part_controllers/connectivity_preserving_model_part_duplication_controller.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/model_part_controllers/connectivity_preserving_model_part_duplication_controller.py)