---
title: Optimization Problem
keywords: 
tags: [optimization problem, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

The ```OptimizationProblem``` is a data container which is used to pass information between components used in the ```OptimizationAnalysis```. There is only one ```OptimizationProblem``` per ```OptimizationAnalysis```

## Data storage

It has two types of data containers.
1. Static data storage.
2. Historical data storage.
3. Component storage.

### Static data storage

Static data storage can only have one value for the whole optimization analysis. Overwriting is not allowed. If required to overwrite, then one should first delete the exiting data, and then write new data. This uses an instance of [BufferedDict](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/utilities/buffered_dict.py).

### Historical data storage

Historical data storage can have different values for different iterations of the whole optimization analysis. Overwriting is not allowed by default (it can be specified). This uses an instance of [BufferedDict](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/utilities/buffered_dict.py).

### Component storage

```OptimizationProblem``` will also store all the instances of ```ExecutionPolicy```s, ```ResponseFunction```s and ```Control```s. Each of the component can have its own data storage within the ```OptimizationProblem```.

The component specific storage can be accessed by the use of [```ComponentDataView```](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/utilities/component_data_view.py). Following python code snippet illustrates a use case
```python
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

model = Kratos.Model()

echo_level = 1
opt_problem = OptimizationProblem(echo_level)

parameters = Kratos.Parameters("""{
    "name"    : "test",
    "type"    : "independent_analysis_execution_policy",
    "settings": {
        "analysis_type"    : "orchestrators.SequentialOrchestrator",
        "analysis_settings": {
            "orchestrator": {
                "settings": {
                    "stage_checkpoints": false
                }
            },
            "stages":[]
        }
    }
}""")
execution_policy = ExecutionPolicyDecorator(model, parameters, opt_problem)

## add the component to opt problem
opt_problem.AddComponent(execution_policy)

# retrieve the component specific data container from opt_problem
data_view = ComponentDataView(execution_policy, opt_problem)

## set buffer to store data for 2 iterations in memory
data_view.SetDataBuffer(2)

## get the historical BufferedDict
historical_data = data_view.GetBufferedData()
print(historical_data)

## get the static BufferedDict
static_data = data_view.GetUnBufferedData()
print(static_data)
```