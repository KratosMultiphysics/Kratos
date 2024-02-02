---
title: Comb
keywords: 
tags: [combing, expressions]
sidebar: kratos_expressions
summary: 
---

## Introduction
This allows combining multiple expressions in to one expression. **The number of entities on all the combining expressions should have the same.**

Example:
All provided expressions in the list_of_expression.
```
For example, let the list_of_expression contain the following expressions:
exp_list[0] = data{1, 2, 3} with 3 items, and item shape = []
                   -  -  -
exp_list[1] = data{4, 5, -1, 6, 7, -1, 8, 9, -1} with 3 items, and item shape = [3]
                   --------  --------  --------
The resulting expression has item shape [3] with 4 items:
output_exp = [1, 4, 5, -1, 2, 6, 7, -1, 3, 8, 9, -1]
                ---------     --------     --------
```
This creates a lazy expression, hence it has a constant cost complexity irrespective of the data size (The expression won't be evaluated unless @ref Expression::Evaluate is called).

## Use case
Following code snippet illustrates how the above example can be replicated.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")

node_1 = model_part.CreateNewNode(1, 0, 0, 0)
node_2 = model_part.CreateNewNode(2, 1, 0, 0)
node_3 = model_part.CreateNewNode(3, 1, 1, 0)

# setting VELOCITY and PRESSURE of each node
node_1.SetValue(Kratos.VELOCITY, Kratos.Array3([4,5,-1]))
node_2.SetValue(Kratos.VELOCITY, Kratos.Array3([6,7,-1]))
node_3.SetValue(Kratos.VELOCITY, Kratos.Array3([8,9,-1]))
node_1.SetValue(Kratos.PRESSURE, 1)
node_2.SetValue(Kratos.PRESSURE, 2)
node_3.SetValue(Kratos.PRESSURE, 3)

# create expressions
u_exp = Kratos.Expression.NodalExpression(model_part)
p_exp = Kratos.Expression.NodalExpression(model_part)

# read data in to expressions
Kratos.Expression.VariableExpressionIO.Read(u_exp, Kratos.VELOCITY, False)
Kratos.Expression.VariableExpressionIO.Read(p_exp, Kratos.PRESSURE, False)

# now combine
combined_exp = Kratos.Expression.Utils.Comb([p_exp, u_exp])

print(combined_exp.Evaluate())
```
Expected output:
```bash
[[ 1.  4.  5. -1.]
 [ 2.  6.  7. -1.]
 [ 3.  8.  9. -1.]]
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
