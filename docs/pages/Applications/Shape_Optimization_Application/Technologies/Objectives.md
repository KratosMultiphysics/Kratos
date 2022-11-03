---
title: Objectives
keywords: 
tags: [Objectives.md]
sidebar: shape_optimization_application
summary: 
---

Objectives are used to determine the shape update in the shape optimization algorithm. The relevant derivatives of objectives are computed in order to calculate the appropriate shape update which will result in an improved objective. Objectives can be defined in the json settings of the problem as shown below.

```json
    "objectives": [
        {
            "identifier": "objective_identifier",
            "type": "minimization",
            "analyzer": "kratos",
            "response_settings": {},
            "project_gradient_on_surface_normals": true
        }
    ],
```

| Option | Allowed values |
| ------------- | ------------- |
| identifier  | Identifier of the objective. This can be a user defined name without spaces.  |
| type  | Type of the objective. It can be either "minimization" or "maximization" [See](#types-of-the-objective). |
| analyzer  | Type of the analyzer. "kratos" or "external" [See](#analyzer).|
| response_settings  | Settings for the response function used to calculate relavant quantities of the objective. [See](#response-function-settings) |
| project_gradient_on_surface_normals  | "true" or "false". Determines whether calculated gradients of objective is projected onto surface normals.|

## Types of the objective

This determines whether the improvement of the objective is expected to be decreasing or increasing. In the case where improvement is expected to be decreasing, then following formulation is used where $$\underline{u} \in \mathbb{R}^M$$ represent state variables obtained by solving governing equations and $$\underline{s} \in \mathbb{R}^N$$ represent design variables which are updated based on the gradient information from the objectives.
<p align="center">$$ \min_{\underline{s} \in\ \mathbb{R}^N} f\left(\underline{u}, \underline{s}\right) $$</p>

In the case where improvement is expected to be increasing then the problem is transformed in to a minimization problem by negating the objective as depicted in the following equation.
<p align="center">$$ \max_{\underline{s} \in\ \mathbb{R}^N} f\left(\underline{u}, \underline{s}\right) =  \min_{\underline{s} \in\ \mathbb{R}^N} -f\left(\underline{u}, \underline{s}\right)$$</p>

Therefore, the standardized value of objective is computed as follows ([source](https://github.com/KratosMultiphysics/Kratos/blob/0048ec0790af5b356039ee4829d78ff0deb2d640/applications/ShapeOptimizationApplication/python_scripts/communicator_factory.py#L222)):
<p align="center">$$ f_{std}\left(\underline{u}, \underline{s}\right) = \begin{cases} -f_{std}\left(\underline{u}, \underline{s}\right) \quad &\textit{if maximization} \\f_{std}\left(\underline{u}, \underline{s}\right) \quad &\textit{if minimization} \end{cases}$$</p>

If the objective is dependent on the state variables (i.e. $$ \underline{u} $$), then the governing equations of the primal problem is imposed as constraints on the optimization problem by forming the Lagrange function as given below where $$\underline{R} = 0$$ are the residuals of the primal giverning equation and $$\lambda$$ are the adjoint variables.

<p align="center">$$ L = f_{std}\left(\underline{u}, \underline{s}\right) + \left(\underline{\lambda}^T\underline{R}\right) $$</p>

Then the total derivative of the objective is taken from the following equation ([source](https://github.com/KratosMultiphysics/Kratos/blob/0048ec0790af5b356039ee4829d78ff0deb2d640/applications/ShapeOptimizationApplication/python_scripts/communicator_factory.py#L239)):
<p align="center">$$ \left(\frac{df}{d\underline{s}}\right)_{std} = \frac{dL}{d\underline{s}} = \frac{\partial L}{\partial \underline{s}} $$</p>

These Lagrange multipliers (i.e. $$ \lambda $$) are computed using the [adjoint approach](../General/Sensitivity_Analysis/Adjoint_approach.html) either fully analytic methodology or semi-analytic methodology.

## Analyzer

Analyzer is used to analyze the given optimization problem and compute respective objective values and gradients to obtain improved objectives as desired by the user. If objectives defined in Kratos is used, then it will invoke kratos internal analyzer. If not user is allowed to use "external" keyword to indicate user defined response functions with user defined analyzer is being used.

## Response function settings

This includes settinsg required to build the response function. Even though these settings may be different for each response function, the first entry "response_type" should be provided with the name of the response function the user expects to use in the respective objective function.

Hello