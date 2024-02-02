---
title: C Array Expression IO
keywords: 
tags: [c array expression io, numpy, scipy, Vector]
sidebar: kratos_expressions
summary: 
---

## Introduction

```CArrayExpressionIO``` is used to exchange data between `C` type arrays and expressions. In python, it can be used to exchange data between numpy arrays, ```Kratos::Vectors```, or any other contiguous array of supported types.

## Reading or viewing a numpy array

Please refer to [Working with numpy](../General/Working_with_Numpy.html) section.

## Writing to numpy array
In MPI runs, expressions will only write data from the ```LocalMesh``` to the numpy array. The size of the numpy array should be compatible with the expression shape. The following code snippet shows an example:
```python
import numpy
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
node_1.SetValue(Kratos.VELOCITY, Kratos.Array3([1,2,3]))
node_2.SetValue(Kratos.VELOCITY, Kratos.Array3([3,4,5]))

# now create the expression:
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# now read the VELOCITY from the non-historical container
Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.VELOCITY, False)

# now get the numpy array from the expression
numpy_array = numpy.zeros((2,3), dtype=numpy.float64)
Kratos.Expression.CArrayExpressionIO.Write(nodal_expression, numpy_array)

# first dimension of the numpy array always represents number of entities in the expression (local mesh entities only)
# rest of the dimensions represent the dimensions of the entity data. In this case, VELOCITY have only three components,
# it shows 3 as the last dimension.
print("numpy array = ", numpy_array)
```

## Reading and writing to Kratos::Vector
```CArrayExpressionIO``` also allows reading and writing from ```Kratos::Vector``` data containers. In this case, the size of the ```Kratos::Vector``` should be the flattened size of the shape of interest which the user will be working on.

The arguments for reading from a ```Kratos::Vector``` are, first the expression should be passed, secondly the ```Kratos::Vector``` should be passed. Finally the shape of the entity data should be passed. Following code snippet illustrates the usage.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)

vector = Kratos.Vector([1,2,3,4,5,6])

# now create the expression:
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# now read in the Kratos::Vector
Kratos.Expression.CArrayExpressionIO.Read(nodal_expression, vector, [3])

# now write the VELOCITY to the non-historical container
Kratos.Expression.VariableExpressionIO.Write(nodal_expression, Kratos.VELOCITY, False)

for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    print(f"node_id: {node.Id}, velocity: [{velocity[0]}, {velocity[1]}, {velocity[2]}]")
```

Writing to a ```Kratos::Vector``` is illustarted by the following code snippet. If the passed vector (to where the expression values are written) is not with the correct size, it will be resized. This will also only write the data values of the local entities in the LocalMesh.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
node_1.SetValue(Kratos.VELOCITY, Kratos.Array3([1,2,3]))
node_2.SetValue(Kratos.VELOCITY, Kratos.Array3([3,4,5]))

# now create the expression:
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# now read the VELOCITY from the non-historical container
Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.VELOCITY, False)

# now write the VELOCITY to the non-historical container
vector = Kratos.Vector()
Kratos.Expression.CArrayExpressionIO.Write(nodal_expression, vector)
print(vector)
```

