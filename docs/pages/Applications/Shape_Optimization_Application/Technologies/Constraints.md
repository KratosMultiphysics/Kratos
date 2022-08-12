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

 Constraints are applied to the optimization problem using lagrange multipliers in the optimization process as depicted in following equation.

Assume following is the objective which is used in the optimization problem.

<p align="center">$$ \min_{\underline{s} \in\ \mathbb{R}^N} f\left(\underline{u}, \underline{s}\right) $$</p>

And the constraints are applied by formulating a Lagrange function including objectives and constraints as depicted below where $$\lambda$$ represents Lagrange multiplier, $$H$$ represents constraint, $$H_0$$ represents the target value of the respective constraint and $$w$$ is the scaling factor.

<p align="center">$$ L = f\left(\underline{u}, \underline{s}\right) + w\lambda\left(H - H_0\right) $$</p>

## Types of constraints

There are mainly five types of constraints.

### Equality constraints

Equality constratints are always active. Therefore, the Lagrange multipliers are always computed for all the equality constraints for all the design iterations in the optimization process. This is specified by specifying "=" as "type" in the constraint settings. Then the value of the constraint which it is equal to is specified by either giving "reference" as "initial_value" or "reference" as "specified_value". In the case of "reference" = "specified_value" then user has to specify "reference_value" field to desired constraint vlaue. Example equality constraint settings block is shown below.

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

Then the value of the constraint which it is compared against to is specified by either giving “reference” as “initial_value” or “reference” as “specified_value”. In the case of “reference” = “specified_value” then user has to specify “reference_value” field to desired constraint vlaue. Example equality constraint settings block is shown below.

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