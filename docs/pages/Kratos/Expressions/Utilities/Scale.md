---
title: Scale
keywords:
tags: [scale, expression]
sidebar: kratos_expressions
summary:
---

## Introduction

This scales value of each component to the specified power given by either a floating value or as an expression.
1. If it is scaled given by a scalar, then all the components are raised to that power.
2. If it is scaled given by a expression and the expression has the same shape as the input expression, then each component is raised to power by the expressions same component given by the second expression.
3. If it is scaled given by a expression and the expression has shape of a scalar expression, then each component is raised to power by the second expression's entity value.

Assume the input expression is given by $$\underline{\mathbb{u}} = \left\lbrace u_{ij},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$ where $$i^{th}$$ entity's $$j^{th}$$ component is represented by $$u_{ij}$$ with $$i\in \left[0, M\right)$$ for each entity and $$j\in \left[0, N\right)$$ for each component in each entity. Following equation illustrates the formulation of the resulting expression.

Case 1:
<p align="center">$$ Scale\left(\underline{\mathbb{u}}, P\right) = \left\lbrace u_{ij}\times P,  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$</p>

Case 2:

The expression with scaling is illustrated by $$\underline{\mathbb{P}} = \left\lbrace p_{ij},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$
<p align="center">$$ Scale\left(\underline{\mathbb{u}}, \underline{\mathbb{P}}\right) = \left\lbrace u_{ij}\times{p_{ij}},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$</p>

Case 3:

The expression with scaling is illustrated by $$\underline{\mathbb{P}} = \left\lbrace p_{i},  \forall i\in\left[0, M\right)\right\rbrace$$
<p align="center">$$ Scale\left(\underline{\mathbb{u}}, \underline{\mathbb{P}}\right) = \left\lbrace u_{ij}\times{p_{i}},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$</p>

## Use cases
Following code snippet illustrates how to use ```Scale```.
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