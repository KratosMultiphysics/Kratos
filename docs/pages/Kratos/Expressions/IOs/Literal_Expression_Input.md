---
title: Literal Expression Input
keywords:
tags: [literal values, expressions, scalars]
sidebar: kratos_expressions
summary:
---

## Introduction

```LiteralExpressionIO``` is used to represent literal values. Following literal types are supported:
* int
* double
* array_1d<double, 3>
* array_1d<double, 4>
* array_1d<double, 6>
* array_1d<double, 9>
* Vector
* Matrix

This IO will create an expression with the specified literal for all entities. Since it is going to be constant for all entities, the literal will be stored only once for
all entities making the memory occupation the least. When the expression is evaluated, this stored literal will be read by all entities.

## Setting literal expressions

Following code snippet shows how to set literal expressions
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)

# create the expression
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# create the literal expression
Kratos.Expression.LiteralExpressionIO.SetData(nodal_expression, Kratos.Array3([1, 2, 3]))

# now write the value to non-historical velocity
Kratos.Expression.VariableExpressionIO.Write(nodal_expression, Kratos.VELOCITY, False)

# now check
for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    print(f"node_id: {node.Id}, velocity = [{velocity[0]}, {velocity[1]}, {velocity[2]}]")
```

## Setting literal expression to zero

The ```LiteralExpressionIO``` also can be used to set entity data to zero for a given variable. **In the case of dynamic data types, it will create empty data sets [Such as Vector with size 0, Matrix with shape (0,0)**. The passed variable is only used to determine the entity data shape. Following code-snippet shows an example:
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)

node_1.SetValue(Kratos.PRESSURE, 10)
node_2.SetValue(Kratos.PRESSURE, 20)

# create the expression
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# read the PRESSURE
Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.PRESSURE, False)

print("PRESSURE before resetting:\n", nodal_expression.Evaluate())

# create the literal expression
Kratos.Expression.LiteralExpressionIO.SetDataToZero(nodal_expression, Kratos.PRESSURE)

# now write the value to non-historical velocity
Kratos.Expression.VariableExpressionIO.Write(nodal_expression, Kratos.TEMPERATURE, False)

# now check
for node in model_part.Nodes:
    print(f"node_id: {node.Id}, pressure = {node.GetValue(Kratos.PRESSURE)}, temperature = {node.GetValue(Kratos.TEMPERATURE)}")
```

## Use cases

Lets assume you want to get inverse of a scalar expressions such as PRESSURE. In that case you can easily use ```LiteralExpressionIO``` to compute it.
Following is an example.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)

node_1.SetValue(Kratos.PRESSURE, 10)
node_2.SetValue(Kratos.PRESSURE, 20)

# create the expression
nodal_expression = Kratos.Expression.NodalExpression(model_part)

# read the PRESSURE
Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.PRESSURE, False)

print("PRESSURE before resetting:\n", nodal_expression.Evaluate())

scalar_expression = Kratos.Expression.NodalExpression(model_part)
Kratos.Expression.LiteralExpressionIO.SetData(scalar_expression, 1.0)
inverse_expression = scalar_expression / nodal_expression

# now write the value to non-historical velocity
Kratos.Expression.VariableExpressionIO.Write(inverse_expression, Kratos.TEMPERATURE, False)

# now check
for node in model_part.Nodes:
    velocity = node.GetValue(Kratos.VELOCITY)
    print(f"node_id: {node.Id}, pressure = {node.GetValue(Kratos.PRESSURE)}, temperature = {node.GetValue(Kratos.TEMPERATURE)}")
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

# now create the expression by reading non historical velocity:
exp = Kratos.Expression.LiteralExpressionIO.Input(model_part, Kratos.Array3([1,2,3]), Kratos.Globals.DataLocation.NodeNonHistorical).Execute()

# do some arithmetic operations
exp *= 2.0

# now write the expression value to model part as ACCELERATION in the historical container
Kratos.Expression.VariableExpressionIO.Output(model_part, Kratos.ACCELERATION, Kratos.Globals.DataLocation.NodeNonHistorical).Execute(exp)

# now print and see
for node in model_part.Nodes:
    acceleration = node.GetValue(Kratos.ACCELERATION)
    print(f"node_id: {node.Id}, acceleration: [{acceleration[0]},{acceleration[1]},{acceleration[2]}]")
```
