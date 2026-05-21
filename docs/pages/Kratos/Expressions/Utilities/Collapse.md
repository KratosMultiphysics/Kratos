---
title: Collapse
keywords: 
tags: [collapse, lazy expression, expressions]
sidebar: kratos_expressions
summary: 
---

## Introduction
```Collapse``` method is used to collapse the lazy expression tree. **This is useful when the lazy expression tree is larger which may occupy or exceeds the RAM
limitations.**.

The resulting expression retains the original's shape, having the same evaluated values. But the resulting expression will not be having the tree structure since it is collapsed.

## Use case
The following example shows how to compute the maximum depth of the expression hierarchy, and collapse it down to a single array.
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

# Destroy the original expression to free up memory
del nodal_expression

print(collapsed_exp.GetMaxDepth())
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
Process Id: 499448 
3
1
```
