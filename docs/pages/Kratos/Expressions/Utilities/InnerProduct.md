---
title: InnerProduct
keywords: 
tags: [inner product, expressions]
sidebar: kratos_expressions
summary: 
---

## Introduction
`InnerProduct `computes component wise inner product between two expressions. Assume the input expression is given by $$\underline{\mathbf{u}} = \left\lbrace u_{ij} |  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$ where the $$i^{th}$$ entity's $$j^{th}$$ component is represented by $$u_{ij}$$ with $$i\in \left[0, M\right)$$ for each entity and $$j\in \left[0, N\right)$$ for each component in each entity. The second input expression is given by $$\underline{\mathbf{v}} = \left\lbrace v_{ij} |  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$ where the $$i^{th}$$ entity's $$j^{th}$$ component is represented by $$u_{ij}$$ with $$i\in \left[0, M\right)$$ for each entity and $$j\in \left[0, N\right)$$ for each component in each entity.

**Both expressions should have the same shape.**

<p align="center">$$ InnerProduct\left(\underline{\mathbf{u}}, \underline{\mathbf{v}}\right) = \sum_{(i,j)\in\left[0, M\right)\times\left[0,N\right)}u_{ij} \times v_{ij}$$</p>

## Use cases
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

print(Kratos.Expression.Utils.InnerProduct(nodal_expression, nodal_expression))
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
Process Id: 513031
285.0
```
