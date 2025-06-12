---
title: Mass Response Function
keywords: 
tags: [mass, response function, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

This response function computes the mass of the given list of model parts.

## Formulation

Following formulation is used.
<p align="center">$$ J   = \sum_{\forall \Omega_i \in \Omega} m_i $$</p>

## Json settings
Following code snippet illustrates json settings used in this response function.
```json
{
    "name": "mass",
    "type": "mass_response_function",
    "settings": {
        "evaluated_model_part_names": [
            "Structure"
        ]
    }
}
```

| Option | Allowed values |
| ------------- | ------------- |
| name | A unique string |
| type  | "mass_response_function"  |
| evaluated_model_part_names  | List of model part names to compute mass |

## Source files

* [applications/OptimizationApplication/python_scripts/responses/mass_response_function.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/responses/mass_response_function.py)
* [Doxygen](doxygen) TODO