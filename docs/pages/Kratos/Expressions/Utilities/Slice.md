---
title: Slice
keywords: 
tags: [slice, expressions]
sidebar: kratos_expressions
summary: 
---

## Introduction
`Slice` method slices existing expression to expressions having a new shape with some sub components from the original expression.
```
Assume an @ref Expression of shape [3] and 3 entities with
the following data in the flattened representation:
@code
data = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        <- 1 ->  <- 2 ->  <- 3 ->
Data for entity 1 is represented with <-1->.

Let Offset = 1 and Stride = 2. The resulting sliced expression
then represents the following data:

output_data = [2, 3, 5, 6, 8, 9]
output container shape = [2] = equal to Stride.
```

## Use case
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

# read in the expression
Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.VELOCITY, False)

# slice
sliced_exp = Kratos.Expression.Utils.Slice(nodal_expression, 1, 2)

print(sliced_exp.Evaluate())
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
Process Id: 540110 
[[2. 3.]
 [5. 6.]
 [8. 9.]]
```
