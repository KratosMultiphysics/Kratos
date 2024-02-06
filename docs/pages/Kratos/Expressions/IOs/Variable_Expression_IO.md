---
title: Variable Expression IO
keywords: 
tags: [variable expression io, expressions, variable]
sidebar: kratos_expressions
summary: 
---

## Introduction

```VariableExpressionIO``` is used to read ```Kratos::Variable``` from a container to an expression vice versa. The expression will only read or write to the entities (i.e. nodes/conditions/elements) in the local mesh. If writing an expression to a variable in a container during an MPI run, synchronization between ranks is automatically done.

### Variable types

The following Kratos variable types are supported:
* ```Kratos::Variable<int>```
* ```Kratos::Variable<double>```
* ```Kratos::Variable<array_1d<double, 3>>```
* ```Kratos::Variable<array_1d<double, 4>>```
* ```Kratos::Variable<array_1d<double, 6>>```
* ```Kratos::Variable<array_1d<double, 9>>```
* ```Kratos::Variable<Vector>```
* ```Kratos::Variable<Matrix>```

### Shape of the resulting expression
If the data type has static size (such as ```Kratos::Variable<double>```, ```Kratos::Variable<array_1d<double, 3>>```), the result's shape is the number of local entities, followed by the static size of the type. However, if the data type has dynamic size (such as ```Kratos::Variable<Vector>```, ```Kratos::Variable<Matrix>```), then all entries in the container are assumed to have the same shape as the first entry. [In MPI runs, the first entry's shape is synchronized across all ranks to determine the final shape of the expression].


## Reading/Writing nodal values

```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
node_1.SetValue(Kratos.VELOCITY, Kratos.Array3([1,2,3]))
node_2.SetValue(Kratos.VELOCITY, Kratos.Array3([3,4,5]))

# now create the expression:
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# now read the VELOCITY from the non-historical container
Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.VELOCITY, False)

# do some arithmetic operations
nodal_expression *= 2.0

# now write the expression value to model part as ACCELERATION in the historical container
Kratos.Expression.VariableExpressionIO.Write(nodal_expression, Kratos.ACCELERATION, True)

# now print and see
for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    acceleration = node.GetSolutionStepValue(Kratos.ACCELERATION)
    print(f"node_id: {node.Id}, velocity: [{velocity[0]},{velocity[1]},{velocity[2]}], acceleration: [{acceleration[0]},{acceleration[1]},{acceleration[2]}]")
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
node_id: 1, velocity: [1.0,2.0,3.0], acceleration: [2.0,4.0,6.0]
node_id: 2, velocity: [3.0,4.0,5.0], acceleration: [6.0,8.0,10.0]
```

## Reading/Writing condition values
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
model_part.CreateNewNode(2, 0.0, 1.0, 0.0)

prop = model_part.CreateNewProperties(1)
cond1 = model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], prop)
cond1.SetValue(Kratos.VELOCITY, Kratos.Array3([1,2,3]))

# now create the expression:
condition_expression = Kratos.Expression.ConditionExpression(model_part)

# now read the VELOCITY from the non-historical container
Kratos.Expression.VariableExpressionIO.Read(condition_expression, Kratos.VELOCITY)

# do some arithmetic operations
condition_expression *= 2.0

# now write the expression value to model part as ACCELERATION in the historical container
Kratos.Expression.VariableExpressionIO.Write(condition_expression, Kratos.ACCELERATION)

# now print and see
for condition in model_part.Conditions:
    velocity = condition.GetValue(Kratos.VELOCITY)
    acceleration = condition.GetValue(Kratos.ACCELERATION)
    print(f"condition_id: {condition.Id}, velocity: [{velocity[0]},{velocity[1]},{velocity[2]}], acceleration: [{acceleration[0]},{acceleration[1]},{acceleration[2]}]")
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
condition_id: 1, velocity: [1.0,2.0,3.0], acceleration: [2.0,4.0,6.0]
```

## Reading/Writing element values
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
model_part.CreateNewNode(3, 1.0, 1.0, 0.0)

prop = model_part.CreateNewProperties(1)
cond1 = model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
cond1.SetValue(Kratos.VELOCITY, Kratos.Array3([1,2,3]))

# now create the expression:
element_expression = Kratos.Expression.ElementExpression(model_part)

# now read the VELOCITY from the non-historical container
Kratos.Expression.VariableExpressionIO.Read(element_expression, Kratos.VELOCITY)

# do some arithmetic operations
element_expression *= 2.0

# now write the expression value to model part as ACCELERATION in the historical container
Kratos.Expression.VariableExpressionIO.Write(element_expression, Kratos.ACCELERATION)

# now print and see
for element in model_part.Elements:
    velocity = element.GetValue(Kratos.VELOCITY)
    acceleration = element.GetValue(Kratos.ACCELERATION)
    print(f"element_id: {element.Id}, velocity: [{velocity[0]},{velocity[1]},{velocity[2]}], acceleration: [{acceleration[0]},{acceleration[1]},{acceleration[2]}]")
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
element_id: 1, velocity: [1.0,2.0,3.0], acceleration: [2.0,4.0,6.0]
```

## Using expressions without the model parts
The ```NodalExpression```, ```ConditionExpression``` and ```ElementExpression``` has an expression which can be directly used if required. The advantage of working
with the ```Expression``` directely is, then it is not bound to a model part of a ```DataValueContainer```. Hence, these expressions can be interchanged if required in
advanced use cases. Following code snippet shows how to use bare ```Expressions```.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
node_1.SetValue(Kratos.VELOCITY, Kratos.Array3([1,2,3]))
node_2.SetValue(Kratos.VELOCITY, Kratos.Array3([3,4,5]))

# now create the expression by reading non historical velocity:
exp = Kratos.Expression.VariableExpressionIO.Input(model_part, Kratos.VELOCITY, Kratos.Globals.DataLocation.NodeNonHistorical).Execute()

# do some arithmetic operations
exp *= 2.0

# now write the expression value to model part as ACCELERATION in the historical container
Kratos.Expression.VariableExpressionIO.Output(model_part, Kratos.ACCELERATION, Kratos.Globals.DataLocation.NodeNonHistorical).Execute(exp)

# now print and see
for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    acceleration = node.GetValue(Kratos.ACCELERATION)
    print(f"node_id: {node.Id}, velocity: [{velocity[0]},{velocity[1]},{velocity[2]}], acceleration: [{acceleration[0]},{acceleration[1]},{acceleration[2]}]")
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
node_id: 1, velocity: [1.0,2.0,3.0], acceleration: [2.0,4.0,6.0]
node_id: 2, velocity: [3.0,4.0,5.0], acceleration: [6.0,8.0,10.0]
```
