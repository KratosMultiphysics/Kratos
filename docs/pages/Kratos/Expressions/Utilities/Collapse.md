---
title: Collapse
keywords: 
tags: [collapse, lazy expression, expressions]
sidebar: kratos_expressions
summary: 
---

## Introduction
```Collapse``` method is used to collapse the lazy expression tree. **This is useful when the lazy expression tree is larger which occupies or exceeds the RAM
limitations.**.

Using this method on an expression which will always result in an expression having the same shape, having the same evaluated values. But the resulting expression will be not having the tree structure since it is collapsed.

## Use case
The max depth of an expression can be checked to identify the max depth of the lazy expression as explained in the following code snippet. It also explains how to collapse the tree hierarchy of the lazy expressions.
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

nodal_expression *= 2
nodal_expression += 1

print(nodal_expression.GetMaxDepth())

collapsed_exp = Kratos.Expression.Utils.Collapse(nodal_expression)
print(collapsed_exp.GetMaxDepth())
```