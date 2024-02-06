---
title: Abs
keywords: 
tags: [abs, expression]
sidebar: kratos_expressions
summary: 
---

## Introduction

This computes the componentwise absolute value of the given expression. Assume the input expression is given by $$\underline{\mathbf{u}} = \left\lbrace u_{ij},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$ where the $$i^{th}$$ entity's $$j^{th}$$ component is represented by $$u_{ij}$$ with $$i\in \left[0, M\right)$$ for each entity and $$j\in \left[0, N\right)$$ for each component in each entity. The Following equation illustrates the formulation of the resulting expression.

<p align="center">$$ Abs\left(\underline{\mathbf{u}}\right) = \left\lbrace \left|u_{ij}\right|,  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$</p>

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
node_1.SetValue(Kratos.VELOCITY, Kratos.Array3([-1,-2,-3]))
node_2.SetValue(Kratos.VELOCITY, Kratos.Array3([-4,-5,-6]))
node_3.SetValue(Kratos.VELOCITY, Kratos.Array3([-7,-8,-9]))

# create the nodal expression
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# read the VELOCITY from non-historical nodal container
Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.VELOCITY, False)

abs_nodal_expression = Kratos.Expression.Utils.Abs(nodal_expression)

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
Process Id: 494037
node_id: 1, velocity=[-1.0, -2.0, -3.0], acceleration = [1.0, 2.0, 3.0]
node_id: 2, velocity=[-4.0, -5.0, -6.0], acceleration = [4.0, 5.0, 6.0]
node_id: 3, velocity=[-7.0, -8.0, -9.0], acceleration = [7.0, 8.0, 9.0]
```

### Using to compute entity wise inf-norm
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

norm_inf_nodal_expression = Kratos.Expression.Utils.EntityMax(Kratos.Expression.Utils.Abs(nodal_expression))

# now write the entity wise norm inf scalar value for checking to the PRESSURE
Kratos.Expression.VariableExpressionIO.Write(norm_inf_nodal_expression, Kratos.PRESSURE, False)

# now printing
for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    pressure = node.GetValue(Kratos.PRESSURE)
    print(f"node_id: {node.Id}, velocity=[{velocity[0]}, {velocity[1]}, {velocity[2]}], pressure = {pressure}")
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
Process Id: 494090
node_id: 1, velocity=[-1.0, -2.0, -3.0], pressure = 3.0
node_id: 2, velocity=[-4.0, -5.0, -6.0], pressure = 6.0
node_id: 3, velocity=[-7.0, -8.0, -9.0], pressure = 9.0
```
