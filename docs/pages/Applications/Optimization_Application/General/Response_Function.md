---
title: Response Function
keywords: 
tags: [response function, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

Response function is used to compute the objective/constraint values and their respective gradients. It purely lives within the physical space, hence all the computed values will be in the physical space.

## Notes
* It will have methods to compute gradients for variables as requested. If a variable is requested, and the gradient computation is not implemented, then it should throw an error.
* ResponseFunction should implement all the variables it depends on. If something is hard to be implemented at the moment of the ```ResponseFunction``` development, then put that variable and throw and error saying not yet implemented.
* ```ResponseFunction``` may or may not use ```ExecutionPolicies```. If it required to have an primal analysis, it must use ```ExecutionPolicy``` to obtain the primal analysis. In this case ```ResponseFunction```. In this case, the ```ResponseFunction::GetEvaluatedModelPart``` will be the domain on which the response is evaluated on and ```ResponseFunction::GetAnalysisModelPart``` will be the domain on which the adjoints may be computed [If no adjoints are used, this method returns None].

## Source files
* [applications/OptimizationApplication/python_scripts/responses/response_function.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/responses/response_function.py)
* [Doxygen](TODO) TODO

