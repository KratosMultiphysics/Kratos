---
title: Expression Arithmetics
keywords: 
tags: [expression, arithmetaics]
sidebar: kratos_expressions
summary: 
---

## Introduction

The ```Expressions``` of type ```NodalExpression```, ```ConditionExpression```, ```ElementExpression``` and their underlying ```Expression``` can be used in arithmetic formulations.

Following operations are supported:
* Addition ```+``` or ```+=```
* Substraction ```-``` or ```-=```
* Multiplication ```*``` or ```*=```
* Division ```/``` or ```/=```
* Power ```**``` or ```**=```

**All these arithmatics always has ```O(1)``` computational cost because, they are not evaluated at the point these arithmetics are used. They are only evaluated when the final value of the ```Expression``` is required.**

## Example usage
Following code snippet explains few use cases for ```NodalExpression```. All the other ```Expression``` also have a similar usage.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0, 0, 0)
node_2 = model_part.CreateNewNode(2, 1, 0, 0)
node_3 = model_part.CreateNewNode(3, 1, 1, 0)

# setting VELOCITY of each node
node_1.SetValue(Kratos.VELOCITY, Kratos.Array3([1,2,3]))
node_2.SetValue(Kratos.VELOCITY, Kratos.Array3([4,5,6]))
node_3.SetValue(Kratos.VELOCITY, Kratos.Array3([7,8,9]))


# create nodal expression
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# read non-historical nodal VELOCITY to expression
Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.VELOCITY, False)

# now we can do arithmetics.
nodal_expression += 1 # Adds 1 to all the components of VELOCITY in all the entities
nodal_expression *= 2 # Multiplies by 2, the all the components of VELOCITY in all the entities
nodal_expression /= 3 # Divides by 2, the all the components of VELOCITY in all the entities
nodal_expression -= 0.5 # Substract by 2, the all the components of VELOCITY in all the entities


# writing the result to the ACCELERATION variable
Kratos.Expression.VariableExpressionIO.Write(nodal_expression, Kratos.ACCELERATION, False)

for node in model_part.Nodes:
    print(node.GetValue(Kratos.ACCELERATION))
```
