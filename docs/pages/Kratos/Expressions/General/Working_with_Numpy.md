---
title: Working with Numpy
keywords: 
tags: [numpy, scipy, expressions, variable expression io, carray expression io]
sidebar: kratos_expressions
summary: 
---

## Introduction

Expressions make working with numpy/scipy or any other thrid party library easier. It allows reading numpy arrays to expressions, and then assign them to Kratos containers using variables.

**When running with MPI, these expressions will hold only data of corresponding containers of the local mesh**.

## Writing to numpy arrays

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

# now get the numpy array from the expression
numpy_nodal_expression = nodal_expression.Evaluate()

# first dimension of the numpy array always represents number of entities in the expression (local mesh entities only)
# rest of the dimensions represent the dimensions of the entity data. In this case, VELOCITY have only three components,
# it shows 3 as the last dimension.
print("Shape of the numpy array = ", numpy_nodal_expression.shape)
```

Expected output:
```console
 |  /           |
 ' /   __| _` | __|  _ \   __|
 . \  |   (   | |   (   |\__ \
_|\_\_|  \__,_|\__|\___/ ____/
           Multi-Physics 9.4."3"-docs/expression_documentation-156476ea1c-Release-x86_64
           Compiled for GNU/Linux and Python3.11 with GCC-13.2
Compiled with threading and MPI support.
Maximum number of threads: 30.
Running without MPI.
Shape of the numpy array =  (2, 3)
```

## Reading from numpy arrays
First create the numpy array independent from Kratos if the number of local entities are known (i.e. number of nodes/conditions/elements in the local mesh). The following numpy data types are supported when reading numpy arrays:
1. ```numpy.int32```
2. ```numpy.float64```

The numpy array can then be **copied** to the expression. The following code snippet shows an example:

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

# now read in the data from numpy array to expression
Kratos.Expression.CArrayExpressionIO.Read(nodal_expression, numpy_array)

# now write the expression data to VELOCITY variable
Kratos.Expression.VariableExpressionIO.Write(nodal_expression, Kratos.VELOCITY, False)

for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    print(f"node id: {node.Id}, velocity: [{velocity[0]}, {velocity[1]}, {velocity[2]}]")
```

Expected output:
```console
 |  /           |
 ' /   __| _` | __|  _ \   __|
 . \  |   (   | |   (   |\__ \
_|\_\_|  \__,_|\__|\___/ ____/
           Multi-Physics 9.4."3"-docs/expression_documentation-156476ea1c-Release-x86_64
           Compiled for GNU/Linux and Python3.11 with GCC-13.2
Compiled with threading and MPI support.
Maximum number of threads: 30.
Running without MPI.
node id: 1, velocity: [1.0, 2.0, 3.0]
node id: 2, velocity: [3.0, 4.0, 5.0]
```

## Moving from numpy arrays

This is similar to [reading from numpy arrays](#reading-from-numpy-arraystest) with the difference that this does not create a copy for the expression.
It uses the memory location of the original numpy array within the expression, eliding the copy and its related memory allocation. An important consequence is that **the original numpy array must be kept in memory for entire lifetime of the expression**, otherwise evaluating the expression will result in a segmentation fault.

Following code snippet illustrates how to move numpy arrays to expressions:
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
Kratos.Expression.CArrayExpressionIO.Move(nodal_expression, numpy_array)

# now write the expression data to VELOCITY variable
Kratos.Expression.VariableExpressionIO.Write(nodal_expression, Kratos.VELOCITY, False)

# now it will print the value in the numpy array
for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    print(f"node id: {node.Id}, velocity: [{velocity[0]}, {velocity[1]}, {velocity[2]}]")

# now change some values in the numpy array
numpy_array[1, 1] = 10.0

# now write the expression data to VELOCITY variable to update the nodal data value container.
Kratos.Expression.VariableExpressionIO.Write(nodal_expression, Kratos.VELOCITY, False)
# now it will print the modified value in the numpy array
for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    print(f"node id: {node.Id}, velocity: [{velocity[0]}, {velocity[1]}, {velocity[2]}]")
```
Expected output:
```console
 |  /           |
 ' /   __| _` | __|  _ \   __|
 . \  |   (   | |   (   |\__ \
_|\_\_|  \__,_|\__|\___/ ____/
           Multi-Physics 9.4."3"-docs/expression_documentation-156476ea1c-Release-x86_64
           Compiled for GNU/Linux and Python3.11 with GCC-13.2
Compiled with threading and MPI support.
Maximum number of threads: 30.
Running without MPI.
node id: 1, velocity: [1.0, 2.0, 3.0]
node id: 2, velocity: [3.0, 4.0, 5.0]
node id: 1, velocity: [1.0, 2.0, 3.0]
node id: 2, velocity: [3.0, 10.0, 5.0]
```
