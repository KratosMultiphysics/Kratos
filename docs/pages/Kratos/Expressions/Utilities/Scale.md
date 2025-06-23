---
title: Scale
keywords: 
tags: [scale, expression]
sidebar: kratos_expressions
summary: 
---

## Introduction

`Scale` scales value of each component to the specified value given by either a floating value or as an expression.
1. If it is scaled by a floating value, then all the components are scaled with the same.
2. If it is scaled by an expression and that expression has the same shape as the first expression, then each component of the first expression is scaled by the second expression's corresponding same component value.
3. If it is scaled by an expression and that expression has shape of a scalar expression, then each component of the first expression is scaled by the second expression's entity value.

Assume the first expression is given by $$\underline{\mathbf{u}} = \left\lbrace u_{ij},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$ where the $$i^{th}$$ entity's $$j^{th}$$ component is represented by $$u_{ij}$$ with $$i\in \left[0, M\right)$$ for each entity and $$j\in \left[0, N\right)$$ for each component in each entity. Following equations illustrate the formulations of the resulting expression.

Case 1:
<p align="center">$$ Scale\left(\underline{\mathbf{u}}, P\right) = \left\lbrace u_{ij}\times P,  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$</p>

Case 2:

The expression with scaling (i.e. second expression) is given by $$\underline{\mathbf{P}} = \left\lbrace p_{ij},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$
<p align="center">$$ Scale\left(\underline{\mathbf{u}}, \underline{\mathbf{P}}\right) = \left\lbrace u_{ij}\times{p_{ij}},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$</p>

Case 3:

The expression with scaling (i.e. second expression) is given by $$\underline{\mathbf{P}} = \left\lbrace p_{i},  \forall i\in\left[0, M\right)\right\rbrace$$
<p align="center">$$ Scale\left(\underline{\mathbf{u}}, \underline{\mathbf{P}}\right) = \left\lbrace u_{ij}\times{p_{i}},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$</p>

## Use cases
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

abs_nodal_expression = Kratos.Expression.Utils.Scale(nodal_expression, 3)

# now write the absolute value for checking to the ACCELERATION
Kratos.Expression.VariableExpressionIO.Write(abs_nodal_expression, Kratos.ACCELERATION, False)

# now printing
for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    acceleration = node.GetValue(Kratos.ACCELERATION)
    print(f"node_id: {node.Id}, velocity=[{velocity[0]}, {velocity[1]}, {velocity[2]}], acceleration = [{acceleration[0]}, {acceleration[1]}, {acceleration[2]}]")
```

Expected output:
```console
 |  /           |
 ' /   __| _` | __|  _ \   __|
 . \  |   (   | |   (   |\__ \
_|\_\_|  \__,_|\__|\___/ ____/
           Multi-Physics 9.4."3"-docs/add_python_processing_locally-eb00abccc7-FullDebug-x86_64
           Compiled for GNU/Linux and Python3.11 with GCC-13.2
Compiled with threading and MPI support.
Maximum number of threads: 30.
Running without MPI.
Process Id: 534910
node_id: 1, velocity=[-1.0, -2.0, -3.0], acceleration = [-3.0, -6.0, -9.0]
node_id: 2, velocity=[-4.0, -5.0, -6.0], acceleration = [-12.0, -15.0, -18.0]
node_id: 3, velocity=[-7.0, -8.0, -9.0], acceleration = [-21.0, -24.0, -27.0]
```
