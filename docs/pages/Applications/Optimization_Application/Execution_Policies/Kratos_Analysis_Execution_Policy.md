---
title: Kratos Analysis Execution Policy
keywords: 
tags: [kratos, execution policy, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

```KratosAnalysisExecutionPolicy``` is used to execute a primal analysis which may or not need to be independent from the model parts of the origin. But they both share the same ```Kratos::Model```.

In this execution policy, when ```Execute``` method is called, following variables in the specified model parts are made to zero.
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
        "analysis_settings": {},
        "analysis_output_settings": {
                    "nodal_solution_step_data_variables": [],
                    "nodal_data_value_variables"        : [],
                    "element_data_value_variables"      : [],
                    "condition_data_value_variables"    : []
        }
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
| analysis_output_settings | Output settings for variables in the model parts listed in `model_part_names`. An [OptimizationProblemVtuOutputProcess](../Processes/Optimization_Problem_Vtu_Output_Process.html) needs to be used to write these fields to files. |
| nodal_solution_step_data_variables | List of nodal solution step variable names |
| nodal_data_value_variables | List of nodal non-historical variable names |
| element_data_value_variables | Element data value variable names |
| condition_data_value_variables | Condition data value variable names |

## Source files
* [Doxygen](TODO) TODO
* [applications/OptimizationApplication/python_scripts/execution_policies/kratos_analysis_execution_policy.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/execution_policies/kratos_analysis_execution_policy.py)
