---
title: Abs
keywords: 
tags: [abs, expression]
sidebar: kratos_expressions
summary: 
---

## Introduction

This computes absolute value of each component of the given expression. Assume the input expression is given by $$\underline{\mathbb{u}} = \left\lbrace u_{ij},  \forall (i,j)\in\left[0, N\right)\times\left[0, M\right)\right\rbrace$$ where $$i^{th}$$ entity's $$j^{th}$$ component is represented by $$u_{ij}$$ with $$i\in \left[0, N\right)$$ for each entity and $$j\in \left[0, M\right)$$. Following equation illustrates the formulation of the resulting expression.

<p align="center">$$ Abs\left(\underline{\mathbb{u}}\right) = \left\lbrace \left|u_{ij}\right|,  \forall (i,j)\in\left[0, N\right)\times\left[0, M\right)\right\rbrace$$</p>

## Use cases

### Using to compute absolute values
Following code snippet illustrates how to use ```Abs```.
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

### Using to compute entity wise inf-norm
Following code snippet illustrates how to use ```Abs``` to compute the entity wise inf-norm.
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
