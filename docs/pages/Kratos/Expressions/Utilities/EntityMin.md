---
title: EntityMin
keywords:
tags: [entity min, min, expressions]
sidebar: kratos_expressions
summary:
---
## Introduction

This get the maximum value from each component in each entity. Assume the input expression is given by $$\underline{\mathbb{u}} = \left\lbrace u_{ij},  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$ where $$i^{th}$$ entity's $$j^{th}$$ component is represented by $$u_{ij}$$ with $$i\in \left[0, M\right)$$ for each entity and $$j\in \left[0, N\right)$$ for each component in each entity. Following equation illustrates the formulation of the resulting expression which is always a scalar expression.

<p align="center">$$ EntityMax\left(\underline{\mathbb{u}}\right) = \left\lbrace \min_{j\in\left[0, N\right)} {u_{ij}},  \forall i\in\left[0, M\right)\right\rbrace$$</p>

## Use cases
Following code snippet illustrates how to use ```EntityMax```.
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

min_nodal_expression = Kratos.Expression.Utils.EntityMin(nodal_expression)

# now write it to PRESSURE for checking.
Kratos.Expression.VariableExpressionIO.Write(min_nodal_expression, Kratos.PRESSURE, False)

# now printing
for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    pressure = node.GetValue(Kratos.PRESSURE)
    print(f"node_id: {node.Id}, velocity=[{velocity[0]}, {velocity[1]}, {velocity[2]}], pressure = {pressure}")
```