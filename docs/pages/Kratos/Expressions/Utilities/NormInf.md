---
title: NormInf
keywords: 
tags: [inf norm, expressions]
sidebar: kratos_expressions
summary: 
---
## Introduction
`NormInf` computes component wise inf norm in expressions. Assume the input expression is given by $$\underline{\mathbf{u}} = \left\lbrace u_{ij} |  \forall (i,j)\in\left[0, M\right)\times\left[0, N\right)\right\rbrace$$ where the $$i^{th}$$ entity's $$j^{th}$$ component is represented by $$u_{ij}$$ with $$i\in \left[0, M\right)$$ for each entity and $$j\in \left[0, N\right)$$ for each component in each entity.

<p align="center">$$ NormInf\left(\underline{\mathbf{u}}\right) = \max_{(i,j)\in\left[0, M\right)\times\left[0,N\right)} \left|u_{ij}\right|$$</p>

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

print(Kratos.Expression.Utils.NormInf(nodal_expression))
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
Process Id: 514980
9.0
```

