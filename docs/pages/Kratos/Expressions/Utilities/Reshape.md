---
title: Reshape
keywords: 
tags: [reshape, expressions]
sidebar: kratos_expressions
summary: 
---

## Introduction
`Reshape` is used to reshape any expression to another new shape. **The total number of components in the new shape should be same as the original shape's total number of components.**

## Use case
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

# now reshape
reshaped_exp = Kratos.Expression.Utils.Reshape(combined_exp, [2, 2])

print(reshaped_exp.Evaluate())
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
Process Id: 528742 
[[[ 1.  4.]
  [ 5. -1.]]

 [[ 2.  6.]
  [ 7. -1.]]

 [[ 3.  8.]
  [ 9. -1.]]]
```
