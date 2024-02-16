---
title: Collective Expressions
keywords: 
tags: [collective expressions, expressions, optimization]
sidebar: optimization_application
summary: 
---


## Introduction

```CollectiveExpressions``` as name suggests is a collection of ```NodalExpression``` and/or ```ConditionExpression``` and/or ```ElementExpression```. You can assume the ```CollectiveExpression``` as a flat vector flattening all the ```Expressions``` it hold. It also allows doing all the arithmetic and other operations which are defined for ```Expressions```. Even though they appear as a flat double vector, A ```CollectiveExpression``` can easily return the underlying ```Expression```s if required.

## Adding and removing ```Expression```s and retrieving expressions.

Following code snippet illustrates how to add and remove ```Expression```s to ```CollectiveExpression```
```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

model = Kratos.Model()
model_part = model.CreateModelPart("test")

nodal_expression = Kratos.Expression.NodalExpression(model_part)
cond_expression = Kratos.Expression.ConditionExpression(model_part)

collective_expression = KratosOA.CollectiveExpression()

# adds a nodal expression and condition expression
collective_expression.Add(nodal_expression)
collective_expression.Add(cond_expression)

# retrieves expressions in the order they were added
for exp in collective_expression.GetContainerExpressions():
    print(exp)

# clears all the expressions
collective_expression.Clear()
```

## Source files
* [Doxygen](TODO) TODO
* [applications/OptimizationApplication/tests/test_collective_expressions.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/tests/test_collective_expressions.py)
