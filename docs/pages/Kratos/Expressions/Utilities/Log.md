---
title: Log
keywords:
tags: [log, expression]
sidebar: kratos_expressions
summary:
---

## Introduction

This computes the component-wise natural logarithm of the given expression. Assume the input expression is given by $$\underline{\mathbf{u}} = \left\lbrace u_{ij},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$ where the $$i^{th}$$ entity's $$j^{th}$$ component is represented by $$u_{ij}$$ with $$i\in \left[0, M\right)$$ for each entity and $$j\in \left[0, N\right)$$ for each component in each entity. The Following equation illustrates the formulation of the resulting expression.

> ##### WARNING
>
> This method returns nan for any component which is $$u_{ij} < 0.0$$ and inf for $$u_{ij} = 0.0$$.

<p align="center">$$ Log\left(\underline{\mathbf{u}}\right) = \left\lbrace log\left(u_{ij}\right),  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$</p>

## Use cases

### Using to compute absolute values
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

# create the nodal expression
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# read the VELOCITY from non-historical nodal container
Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.VELOCITY, False)

log_nodal_expression = Kratos.Expression.Utils.Log(nodal_expression)

# now write the absolute value for checking to the ACCELERATION
Kratos.Expression.VariableExpressionIO.Write(log_nodal_expression, Kratos.ACCELERATION, False)

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
           Multi-Physics 9.5."0"-core/expression/feature/add_log_expression-c757d39762-Release-x86_64
           Compiled for GNU/Linux and Python3.12 with GCC-14.1
Compiled with threading and MPI support.
Maximum number of threads: 36.
Running without MPI.
node_id: 1, velocity=[1.0, 2.0, 3.0], acceleration = [0.0, 0.6931471805599453, 1.0986122886681098]
node_id: 2, velocity=[4.0, 5.0, 6.0], acceleration = [1.3862943611198906, 1.6094379124341003, 1.791759469228055]
node_id: 3, velocity=[7.0, 8.0, 9.0], acceleration = [1.9459101490553132, 2.0794415416798357, 2.1972245773362196]
```
