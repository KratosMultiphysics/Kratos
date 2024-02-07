---
title: Visualization
keywords:
tags: [visualization expression vtu vtk]
sidebar: kratos_expressions
summary:
---

## Introduction

Expressions can be written to the `Vtu` format and can be easily visualized using any compatible software such as [Paraview](https://www.paraview.org/).

This output is compatible with **shared memory** and **distributed memory** parallelization. In **shared memory** parallelization, it will write only one file having the extension `*.vtu`. In **distributed memory** parallelization, it will write one `*.vtu` file per rank, and then one `*.pvtu` for all ranks linking all the `*.vtu` files from all ranks.

## Use case
```python
import numpy
# create the numpy expression. This does not need Kratos as long as you know correctly
# the number of local entities (i.e. nodes/conditions/elements) in the model part.
# it is a good practice to specify the "dtype" in here.
# here also, first dimension represent the number of local entities, rest of the
# dimensions represent the dimensionality of the entity data.
numpy_array = numpy.array([
                    [1,2,3],
                    [3,4,5]], dtype=numpy.float64)


import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)

# now create the expression:
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# now move in the data from numpy array to expression (linking nodal expression with numpy array)
Kratos.Expression.CArrayExpressionIO.Read(nodal_expression, numpy_array)

# creates the vtu output. Needs to be same as the model part used to create the expression
vtu_output = Kratos.VtuOutput(model_part)

# add the expression. First argument is the name to be used in the visualization. second argument is the expression
vtu_output.AddContainerExpression("numpy_array", nodal_expression)

# print the output. Argument specifies the file name prefix.
# this will only print the nodes, since the expression is a nodal expression without any elements.
# if the model part has elements or conditions, then it will print the corresponding geometries with nodal values given in expression.
vtu_output.PrintOutput("output")
```
