---
title: EntitySum
keywords:
tags: [entity sum, sum, expressions]
sidebar: kratos_expressions
summary:
---
## Introduction

This get the sum of value from each component in each entity. Assume the input expression is given by $$\underline{\mathbf{u}} = \left\lbrace u_{ij},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$ where $$i^{th}$$ entity's $$j^{th}$$ component is represented by $$u_{ij}$$ with $$i\in \left[0, M\right)$$ for each entity and $$j\in \left[0, N\right)$$ for each component in each entity. Following equation illustrates the formulation of the resulting expression which is always a scalar expression.

<p align="center">$$ EntitySum\left(\underline{\mathbf{u}}\right) = \left\lbrace v_i, \text{where} \quad v_i = \sum_{j\in\left[0, N\right)} {u_{ij}} \quad  \forall i\in\left[0, M\right)\right\rbrace$$</p>

## Use cases
Following code snippet illustrates how to use ```EntitySum```.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")

node_1 = model_part.CreateNewNode(1, 0, 0, 0)
node_2 = model_part.CreateNewNode(2, 1, 0, 0)
node_3 = model_part.CreateNewNode(3, 1, 1, 0)

# setting VELOCITY of each node
node_1.SetValue(Kratos.VELOCITY, Kratos.Array3([-1,-2,-3]))
node_2.SetValue(Kratos.VELOCITY, Kratos.Array3([-4,-5,-6]))
node_3.SetValue(Kratos.VELOCITY, Kratos.Array3([-7,-8,-9]))

# create the nodal expression
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# read the VELOCITY from non-historical nodal container
Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.VELOCITY, False)

min_nodal_expression = Kratos.Expression.Utils.EntitySum(nodal_expression)

# now write it to PRESSURE for checking.
Kratos.Expression.VariableExpressionIO.Write(min_nodal_expression, Kratos.PRESSURE, False)

# now printing
for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    pressure = node.GetValue(Kratos.PRESSURE)
    print(f"node_id: {node.Id}, velocity=[{velocity[0]}, {velocity[1]}, {velocity[2]}], pressure = {pressure}")
```
Expected output:
```bash
node_id: 1, velocity=[-1.0, -2.0, -3.0], pressure = -6.0
node_id: 2, velocity=[-4.0, -5.0, -6.0], pressure = -15.0
node_id: 3, velocity=[-7.0, -8.0, -9.0], pressure = -24.0
 |  /           |
 ' /   __| _` | __|  _ \   __|
 . \  |   (   | |   (   |\__ \
_|\_\_|  \__,_|\__|\___/ ____/
           Multi-Physics 9.4."3"-docs/expression_documentation-6de5f1a499-Release-x86_64
           Compiled for GNU/Linux and Python3.11 with GCC-13.2
Compiled with threading and MPI support.
Maximum number of threads: 30.
Running without MPI.
```
