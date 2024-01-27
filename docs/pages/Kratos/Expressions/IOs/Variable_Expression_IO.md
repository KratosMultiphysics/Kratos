---
title: Variable Expression IO
keywords: 
tags: [variable expression io, expressions, variable]
sidebar: kratos_expressions
summary: 
---

## Introduction

```VariableExpressionIO``` is used to read ```Kratos::Variable``` data in a data value container to a expression or write an expression data in to a ```Kratos::Variable``` in a specified data value container. The expression will only read or write to the entities (i.e. nodes/conditions/elements) in the local mesh. In the case of writing expression data to a variable in data value container, if it is used in ***distributed memory** parallelized environment, proper synchronization will be automatically done. Hence, all of these operations are compatible with **shared memory** and **distributed memory** parallelizations.

### Variable types

Following variable types are allowed to be read in or write to from an expression.
* ```Kratos::Variable<int>```
* ```Kratos::Variable<double>```
* ```Kratos::Variable<array_1d<double, 3>>```
* ```Kratos::Variable<array_1d<double, 4>>```
* ```Kratos::Variable<array_1d<double, 6>>```
* ```Kratos::Variable<array_1d<double, 9>>```
* ```Kratos::Variable<Vector>```
* ```Kratos::Variable<Matrix>```

### Shape of the resulting expression
In the case of static data types given by static data type variables (such as ```Kratos::Variable<double>```, ```Kratos::Variable<array_1d<double, 3>>```), the shape of the expression will be always have number of local entities in the first dimension, rest will be followed by the static size of the data type. If the variable type being used is of dynamic data type (such as ```Kratos::Variable<Vector>```, ```Kratos::Variable<Matrix>```), then it will represent still the number of local entities in the first dimension. Rest of the dimensions will represent the dimenions of the first data value found in the respective container. [In the case of **distributed memory** parallelized runs, first data value shape synchronization is done to identify the shape of the expression].


## Reading/Writing nodal values
Reading and writing to nodal values have the same interface within ```VariableExpressionIO```. The first argument shall be a `NodalExpression`, the second argumen shall be a ```Kratos::Variable``` and the last argument shall be a boolean indicating true if it is required to read from historical nodal data value container or false if it is required read from non-historical nodal data value container. In the case of **distributed memory** parallelized runs, a proper synchronization between ghost nodes are done when writing is executed.

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

## Reading/Writing condition values
Reading and writing to condition values have the same interface within ```VariableExpressionIO```. The first argument shall be a `ConditionExpression` and the second argumen shall be a ```Kratos::Variable```.
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

## Reading/Writing element values
Reading and writing to element values have the same interface within ```VariableExpressionIO```. The first argument shall be a `ElementExpression` and the second argumen shall be a ```Kratos::Variable```.
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
