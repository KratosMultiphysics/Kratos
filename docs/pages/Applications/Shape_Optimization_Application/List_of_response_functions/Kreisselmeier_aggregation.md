---
title: Kreisselmeier aggregation
keywords: 
tags: [Kreisselmeier_aggregation.md]
sidebar: shape_optimization_application
summary: 
---

## Introduction

Kreisselmeier aggregation is objective with an aggregation methodology used when typical objectives and constraints illustrate ocsillating behaviour when or $$max$$ functions are used. This aggregation methodology generates a smooth functional mimicking the $$max$$ function without the ocsillating behaviour which is good for optimization problems when used with their gradient information.

## Formulation

The response function is computed as follows where $$\beta$$ is chosen value for each entity in the $$\Omega$$ ($$\Omega$$ can be nodal value or elemental gauss point value):
<p align="center">$$ J   = \frac{1}{\rho}\log\left({\sum_\Omega{e^{\rho \frac{\beta}{\alpha}}}}\right)$$</p>

The coefficients are listed in the following table:

|Coefficient | Setting | Description |
|------------|---------|-------------|
|$$\rho$$    | `Aggregation penalty` | This works as a smoothning factor |
|$$\alpha$$ | `Scaling factor` | This works as a scaling factor. Too law values will result in infinity and too large values will result in less weighting in the sensitivity computation. It also accepts `max` as a keyword. Then the solver will compute the max $$\beta$$ in the initial geometry and then use it as a scaling factor. |


## Usage in Kratos Multiphysics

```json
    {
        "gauss_point_values_scalar_variable"            : "PLEASE_SPECIFY_A_SCALAR_VARIABLE",
        "gauss_point_gradient_matrix_variable"          : "PLEASE_SPECIFY_A_MATRIX_VARIABLE",
        "gauss_point_first_derivatives_matrix_variable" : "PLEASE_SPECIFY_A_MATRIX_VARIABLE",
        "gauss_point_second_derivatives_matrix_variable": "PLEASE_SPECIFY_A_MATRIX_VARIABLE",
        "gauss_point_shape_derivatives_matrix_variable" : "PLEASE_SPECIFY_A_MATRIX_VARIABLE",
        "gauss_point_value_scaling_factor"              : "max",
        "aggregation_penalty"                           : 50.0,
        "echo_level"                                    : 0,
        "gradient_mode"                                 : "semi_analytic",
        "gradient_mode_settings": {
            "perturbation_variable_name"           : "PLEASE_SPECIFY_PERTURBATION_SCALAR_VARIABLE_NAME",
            "perturbation_size"                    : 1e-6,
            "design_variable_name_storage_variable": "PLEASE_SPECIFY_STRING_VARIABLE"
        }
    }
```

Following table explains the each settings present in the json.

|Setting|Description|
|-------|-----------|
|gauss_point_values_scalar_variable| The scalar variable to compute the $$\beta$$ value|
|gauss_point_gradient_matrix_variable| The matrix variable to compute gradients partial derivative of $$\beta$$ value|
|gauss_point_first_derivatives_matrix_variable| The matrix variable to compute first partial derivatives of the $$\beta$$ value|
|gauss_point_second_derivatives_matrix_variable| The matrix variable to compute second partial derivatives of the $$\beta$$ value|
|gauss_point_shape_derivatives_matrix_variable| The matrix variable to compute shape partial derivatives of the $$\beta$$ value|
|gauss_point_value_scaling_factor| $$\alpha$$ value  for scaling of the $$\beta$$||
|aggregation_penalty| $$\rho$$ value |
|gradient_mode|Only supports "semi_analytic"|

Following table illustrates "gradient_mode_settings".

|Settings|Description|
|--------|-----------|
|perturbation_variable_name| Name of the variable to be perturbed for finite difference method|
|perturbation_size| Perturbation size of the value|
|design_variable_name_storage_variable| Place holder to store the design variable|

## Source

Location: ["applications/ShapeOptimizationApplication/custom_response_functions/gauss_point_kreisselmeier_aggregation_response_function.h"](https://github.com/KratosMultiphysics/Kratos/blob/shapeopt/kreisselmeier_aggregation/applications/ShapeOptimizationApplication/custom_response_functions/gauss_point_kreisselmeier_aggregation_response_function.h)
