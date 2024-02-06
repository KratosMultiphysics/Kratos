---
title: Nodal Position Expression IO
keywords: 
tags: [nodal position, initial, current, expressions]
sidebar: kratos_expressions
summary: 
---

## Introduction

```NodalPositionExpressionIO``` allows writing and reading nodal positions of a given model part. This can be used on following the configurations:
* Initial -> initial positions of the node
* Current -> positions of the nodes in the current configuration

In the case of **distributed memory** parallelized runs, the reading will only read the nodal positions of the nodes in the local mesh, and the writing will write to local mesh nodes and will afterwards synchronize to update the ghost meshes.

## Reading and writing nodal positions
Following code snippet shows how to read initial and current configuration from nodes.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)

delta = 10.0
node_1.X += delta
node_1.Y += delta

# now create the expression:
initial_nodal_expression = Kratos.Expression.NodalExpression(model_part)
current_nodal_expression = Kratos.Expression.NodalExpression(model_part)

# now read the nodal positions
Kratos.Expression.NodalPositionExpressionIO.Read(initial_nodal_expression,Kratos.Configuration.Initial)
Kratos.Expression.NodalPositionExpressionIO.Read(current_nodal_expression,Kratos.Configuration.Current)

print("initial:\n", initial_nodal_expression.Evaluate())
print("current:\n", current_nodal_expression.Evaluate())

# now do some arithmetics with the exps
initial_nodal_expression *= 2
current_nodal_expression *= 2

# now write the modified nodal positions back to nodes
Kratos.Expression.NodalPositionExpressionIO.Write(initial_nodal_expression,Kratos.Configuration.Initial)
Kratos.Expression.NodalPositionExpressionIO.Write(current_nodal_expression,Kratos.Configuration.Current)

# now read nodal positions
for node in model_part.Nodes:
    print(f"node_id: {node.Id}, initial: [{node.X0}, {node.Y0}, {node.Z0}], current: [{node.X}, {node.Y}, {node.Z}]")
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
initial:
 [[0. 0. 0.]
 [0. 1. 0.]]
current:
 [[10. 10.  0.]
 [ 0.  1.  0.]]
node_id: 1, initial: [0.0, 0.0, 0.0], current: [20.0, 20.0, 0.0]
node_id: 2, initial: [0.0, 2.0, 0.0], current: [0.0, 2.0, 0.0]
```

## Using expressions without the model parts
The ```NodalExpression```, ```ConditionExpression``` and ```ElementExpression``` has an expression which can be directly used if required. The advantage of working
with the ```Expression``` directly is, then it is not bound to a model part of a ```DataValueContainer```. Hence, these expressions can be interchanged if required in
advanced use cases. Following code snippet shows how to use bare ```Expressions```.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)

# now create the expression by reading initial configuration:
exp = Kratos.Expression.NodalPositionExpressionIO.Input(model_part, Kratos.Configuration.Initial).Execute()

# do some arithmetic operations
exp *= 2.0

# now write the expression value to current configuration
Kratos.Expression.NodalPositionExpressionIO.Output(model_part, Kratos.Configuration.Current).Execute(exp)

# now read nodal positions
for node in model_part.Nodes:
    print(f"node_id: {node.Id}, initial: [{node.X0}, {node.Y0}, {node.Z0}], current: [{node.X}, {node.Y}, {node.Z}]")
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
node_id: 1, initial: [0.0, 0.0, 0.0], current: [0.0, 0.0, 0.0]
node_id: 2, initial: [0.0, 1.0, 0.0], current: [0.0, 2.0, 0.0]
```

