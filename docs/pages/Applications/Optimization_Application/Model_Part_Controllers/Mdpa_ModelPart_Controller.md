---
title: Mdpa ModelPart Controller
keywords: 
tags: [mdpa, model part controller, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

This ```MdpaModelPartController``` is used to read model parts from mdpa files.


## Json settings
Following json-snippet illustrates an example use case
```json
{
    "type": "mdpa_model_part_controller",
    "module": "KratosMultiphysics.OptimizationApplication.model_part_controllers",
    "settings": {
        "model_part_name": "Structure",
        "domain_size": 3,
        "input_filename": "Structure"
    }
}
```

| Option | Allowed values |
| ------------- | ------------- |
| type  | "mdpa_model_part_controller"  |
| module  | "KratosMultiphysics.OptimizationApplication.model_part_controllers"  |
| model_part_name  | Model part name to be used. If not found, it will be created |
| domain_size  | Dimensionality of the mesh. Either 2 or 3 |
| input_filename  | mdpa file name without the .mdpa extension |

## Source files
* [Doxygen](TODO) TODO
* [applications/OptimizationApplication/python_scripts/model_part_controllers/mdpa_model_part_controller.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/model_part_controllers/mdpa_model_part_controller.py)

