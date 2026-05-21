---
title: Independent Analysis Execution Policy
keywords: 
tags: [independent, execution policy, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

```IndependentAnalysisExecutionPolicy``` is used to execute a primal analysis which needs to be independent from the model parts of the origin. But they both share the same ```Kratos::Model```.

## Json settings
Following json-snippet illustrates an example use case
```json
{
    "type": "independent_analysis_execution_policy",
    "module": "KratosMultiphysics.OptimizationApplication.execution_policies",
    "settings": {
        "analysis_module"         : "KratosMultiphysics",
        "analysis_type"           : "",
        "analysis_model_part_name": "",
        "analysis_settings"       : {}
    }
}
```

| Option | Allowed values |
| ------------- | ------------- |
| type  | "independent_analysis_execution_policy"  |
| module  | "KratosMultiphysics.OptimizationApplication.model_part_controllers" |
| analysis_module | Where to find the module for `analysis_type`. |
| analysis_type | Analysis type to used to solve the primal analysis. |
| analysis_settings | Settings for the analysis |

## Source files
* [Doxygen](TODO) TODO
* [applications/OptimizationApplication/python_scripts/execution_policies/independent_analysis_execution_policy.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/execution_policies/independent_analysis_execution_policy.py)