---
title: Stepping Analysis Execution Policy
keywords: 
tags: [stepping kratos, execution policy, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

```SteppingAnalysisExecutionPolicy``` is used to execute a primal analysis which may or not need to be independent from the model parts of the origin. But they both share the same ```Kratos::Model```.

In this execution policy, when ```Execute``` method is called, following variables in the specified model parts are made to zero before and thereafter, original values are restored after running the primal analysis.
* STEP
* TIME
* DELTA_TIME

## Json settings
Following json-snippet illustrates an example use case
```json
{
    "type": "independent_analysis_execution_policy",
    "module": "KratosMultiphysics.OptimizationApplication.execution_policies",
    "settings": {
        "model_part_names" : [],
        "analysis_module"  : "KratosMultiphysics",
        "analysis_type"    : "",
        "analysis_settings": {}
    }
}
```

| Option | Allowed values |
| ------------- | ------------- |
| type  | "independent_analysis_execution_policy"  |
| module  | "KratosMultiphysics.OptimizationApplication.model_part_controllers" |
| model_part_names | Model part names to be used for outputting data and resetting variables. |
| analysis_module | Where to find the module for `analysis_type`. |
| analysis_type | Analysis type to used to solve the primal analysis. |
| analysis_settings | Settings for the analysis |

## Source files
* [Doxygen](TODO) TODO
* [applications/OptimizationApplication/python_scripts/execution_policies/stepping_analysis_execution_policy.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/execution_policies/stepping_analysis_execution_policy.py)