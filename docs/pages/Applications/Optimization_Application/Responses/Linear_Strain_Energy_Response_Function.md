---
title: Linear Strain Energy Response Function
keywords: 
tags: [strain energy, response function, optimization]
sidebar: optimization_application
summary: 
---
## Introduction

This computes the summation of strain energy from each element as the response value and its gradient.

## Formulation

Following formulation is used to compute the summation of strain energy from each element in the chosen model part where $$\underline{u}$$ is the displacement vector for the element and $$\mathbf{K}$$ is the stiffness matrix of the element.
<p align="center">$$ J   = \frac{1}{2}\sum_{\forall \Omega_i \in \Omega} \underline{u}^T_i \mathbf{K}_i \underline{u}_i  $$</p>

## Json settings
Following code snippet illustrates json settings used in this response function.
```json
{
    "name": "strain_energy",
    "type": "linear_strain_energy_response_function",
    "settings": {
        "evaluated_model_part_names": [
            "Structure"
        ],
        "primal_analysis_name": "Structure_static",
        "perturbation_size": 1e-8
    }
}
```

| Option | Allowed values |
| ------------- | ------------- |
| name | A unique string |
| type  | "linear_strain_energy_response_function"  |
| evaluated_model_part_names  | List of model part names to compute strain energy |
| primal_analysis_name | Name of the analysis from the list of analysis to be used for strain energy computation |
| perturbation_size | Perturbation size to be used in the semi-analytic method |

## Source files

* [applications/OptimizationApplication/python_scripts/responses/linear_strain_energy_response_function.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/responses/linear_strain_energy_response_function.py)
* [Doxygen](doxygen) TODO
