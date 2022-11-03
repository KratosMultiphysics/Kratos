---
title: Constraints
keywords: 
tags: [Constraints.md]
sidebar: shape_optimization_application
summary: 
---

Constraints are used to constraint a given optimization problem. These constraints can be linear or non-linear. Depending on whether they are linear or non-linear, they may require to solve for one additional adjoint problem to obtain respective sensitivities.

Following snippet illustrates available options for constraints.
```json
"constraints": [
    {
            "identifier"                          : "NO_IDENTIFIER_SPECIFIED",
            "type"                                : "<",
            "scaling_factor"                      : 1.0,
            "reference"                           : "initial_value",
            "reference_value"                     : 1.0,
            "analyzer"                            : "external",
            "response_settings"                   : {},
            "project_gradient_on_surface_normals" : false
    }
]
```

Following table illustrates descriptions of each settings available in each constraint.

| Option | Allowed values |
| ------------- | ------------- |
| identifier  | Identifier of the constraint. This can be a user defined name without spaces.  |
| type  | Type of the constraint. It can be "=", "<", ">", "<=", ">=" [See](#types-of-the-constraint). |
| analyzer  | Type of the analyzer. "kratos" or "external" [See](#analyzer).|
| response_settings  | Settings for the response function used to calculate relavant quantities of the constraint. [See](#response-function-settings) |
| project_gradient_on_surface_normals  | "true" or "false". Determines whether calculated gradients of constraint is projected onto surface normals.|
|scaling_factor| Scaling factor to be used in Lagraneg formulation|
|reference| The reference value type of the response. It can be either "initial_value", "specified_value"|
|referece_value| In the case of "reference_type"="specified_value", then this value is used as the user specified value|

## Constraint formulation

These constraints can be active or inactive depending of the type of the constraint (see [Types of constraints](#types-of-constraints) which illustrates how they can be specified).

Firstly the constraint value is standardized as given below (i.e. objective value is $$g$$, state variables are $$ \underline{w} $$, design variable is $$ s $$ and reference value for constraint is $$ g_{ref} $$) ([source](https://github.com/KratosMultiphysics/Kratos/blob/0048ec0790af5b356039ee4829d78ff0deb2d640/applications/ShapeOptimizationApplication/python_scripts/communicator_factory.py#L227)):
<p align="center">$$ g_{std}\left(\underline{u}, s\right) = \begin{cases} g_{ref} - g\left(\underline{u}, s\right)  \quad &\textit{if "type" = ">" or ">="}\\ g\left(\underline{u}, s\right) - g_{ref} \quad &\textit{otherwise} \end{cases}$$</p>

Then an in-equality constraint is deemed active if $$ g_{std} > 0.0 $$ otherwise they deemed as inactive. Equality constraints are always deemed as active. If they are deemed active, then following equation is used to compute the standardized sensitivities of the constraint ([source](https://github.com/KratosMultiphysics/Kratos/blob/0048ec0790af5b356039ee4829d78ff0deb2d640/applications/ShapeOptimizationApplication/python_scripts/communicator_factory.py#L239)).
<p align="center">$$ \left(\frac{dg}{d\underline{s}}\right)_{std}\left(\underline{u}, s\right) = \begin{cases} - \frac{dg}{d\underline{s}}\left(\underline{u}, s\right)  \quad &\textit{if "type" = ">" or ">="}\\ \frac{dg}{d\underline{s}}\left(\underline{u}, s\right) \quad &\textit{otherwise} \end{cases}$$</p>

if the $$g$$ is dependent on state variables (i.e. $$\underline{u}$$) then, the residuals of the primal governing equations are applied as additional constraints to the $$g$$ using the Lagrange multipliers (i.e. $$\lambda$$). Then the sensitivities of the constraints are computed by formulating a Lagrange function including constraint and its primal governing equation residuals (i.e. $$R$$) as depicted below.

<p align="center">$$ L = g\left(\underline{u}, \underline{s}\right) + \left(\underline{\lambda}^T\underline{R}\right) $$</p>

Then the total derivative of the constraint is taken from the following equation:
<p align="center">$$ \frac{dg}{d\underline{s}} = \frac{dL}{d\underline{s}} = \frac{\partial L}{\partial \underline{s}} $$</p>

These Lagrange multipliers (i.e. $$ \lambda $$) are computed using the [adjoint approach](../General/Sensitivity_Analysis/Adjoint_approach.html) either fully analytic methodology or semi-analytic methodology.

## Types of constraints

There are mainly five types of constraints.

### Equality constraints

Equality constratints are always active. Therefore, the Lagrange multipliers are always computed for all the equality constraints for all the design iterations in the optimization process. This is specified by specifying "=" as "type" in the constraint settings. Then the value of the constraint which it is equal (i.e. $$g_{ref}$$) to is specified by either giving "reference" as "initial_value" or "reference" as "specified_value". In the case of "reference" = "specified_value" then user has to specify "reference_value" field to desired constraint vlaue. Example equality constraint settings block is shown below.

```json
"constraints": [
    {
            "identifier"                          : "NO_IDENTIFIER_SPECIFIED",
            "type"                                : "=",
            "scaling_factor"                      : 1.0,
            "reference"                           : "initial_value",
            "reference_value"                     : 1.0,
            "analyzer"                            : "external",
            "response_settings"                   : {},
            "project_gradient_on_surface_normals" : false
    }
]
```

### In-equality constraints

In-equality constraints are only active when they are violated. Therefore the Lagrange multipliers corresponding to in-equality constraints are only comptued when they are violated. In-equality constraints can be specified by specifying the "type" of the constraint settings as given in the following table.

|Type|Description|
|"<"| Less than the "initial_value" or "specified_value"|
|">"| Greater than the "initial_value" or "specified_value"|
|"<="| Less than or equal to the "initial_value" or "specified_value"|
|">="| Greater than or equal to the "initial_value" or "specified_value"|

Then the value of the constraint which it is compared against (i.e. $$g_{ref}$$) to is specified by either giving “reference” as “initial_value” or “reference” as “specified_value”. In the case of “reference” = “specified_value” then user has to specify “reference_value” field to desired constraint vlaue. Example equality constraint settings block is shown below.

```json
"constraints": [
    {
            "identifier"                          : "NO_IDENTIFIER_SPECIFIED",
            "type"                                : ">=",
            "scaling_factor"                      : 1.0,
            "reference"                           : "specified_value",
            "reference_value"                     : 1.0,
            "analyzer"                            : "external",
            "response_settings"                   : {},
            "project_gradient_on_surface_normals" : false
    }
]
```