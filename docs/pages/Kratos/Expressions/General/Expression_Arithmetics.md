---
title: Expression Arithmetics
keywords: 
tags: [expression, arithmetics]
sidebar: kratos_expressions
summary: 
---

## Introduction

The ```Expressions``` of type ```NodalExpression```, ```ConditionExpression```, ```ElementExpression``` and their underlying ```Expression``` can be used in arithmetic formulations.

Following operations are supported:
* Addition ```+``` or ```+=```
* Subtraction ```-``` or ```-=```
* Multiplication ```*``` or ```*=```
* Division ```/``` or ```/=```
* Power ```**``` or ```**=```

These operations are carried out only when their results are required ([lazy evaluation](https://en.wikipedia.org/wiki/Lazy_evaluation)).

## Example usage
Here is an example for using ```NodalExpression``` (other expression types can be used in a similar fashion).
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
nodal_expression -= 0.5 # Subtract by 2, the all the components of VELOCITY in all the entities

# writing the result to the ACCELERATION variable
Kratos.Expression.VariableExpressionIO.Write(nodal_expression, Kratos.ACCELERATION, False)

for node in model_part.Nodes:
    print(node.GetValue(Kratos.ACCELERATION))
```
Expected output:
```console
 |  /           |
 ' /   __| _` | __|  _ \   __|
 . \  |   (   | |   (   |\__ \
_|\_\_|  \__,_|\__|\___/ ____/
           Multi-Physics 9.4."3"-docs/expression_documentation-156476ea1c-Release-x86_64
           Compiled for GNU/Linux and Python3.11 with GCC-13.2
Compiled with threading and MPI support.
Maximum number of threads: 30.
Running without MPI.
[3](0.833333,1.5,2.16667)
[3](2.83333,3.5,4.16667)
[3](4.83333,5.5,6.16667)
```
