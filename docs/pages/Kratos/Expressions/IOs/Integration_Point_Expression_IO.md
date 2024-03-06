---
title: Integration Point Expression IO
keywords: 
tags: [Integration_Point_Expression_IO.md]
sidebar: kratos_expressions
summary: 
---

## Introduction

```IntegrationPointExpressionIO``` is used to set or get quantities at integration points from elements or conditions.

### Variable types
Following data types are supported for integration point quantity calculations.
* ```int```
* ```double```
* ```array_1d<double, 3>```
* ```array_1d<double, 4>```
* ```array_1d<double, 6>```
* ```array_1d<double, 9>```
* ```Vector```
* ```Matrix```

### Shape of resulting expressions
The resulting expressions computed by ```IntegrationPointExpressionIO``` will have two additional dimensions on top of the data type's dimensions. The first one indicates the number of entities in the expression, while the second dimension is the number of gauss points computed in each entitiy. The rest of the dimensions are carried over from the integration point quantity.

## Reading and writing integration point data
Following code snippet shows how to read and write integration point data.
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

# now read the VELOCITY from the container
Kratos.Expression.IntegrationPointExpressionIO.Read(element_expression, Kratos.VELOCITY)

# the shape of the element expression will be [1, 1, 3] where first "1" is number of elements, second "1" is number of gauss points.
shape = element_expression.Evaluate().shape

# do some arithmetic operations
element_expression *= 2.0

# now read the ACCELERATION from the container
Kratos.Expression.IntegrationPointExpressionIO.Write(element_expression, Kratos.ACCELERATION)
```

## Using expressions without the model parts
The underlying `Expression` instance of ```NodalExpression```, ```ConditionExpression``` and ```ElementExpression``` can be accessed directly, if necessary. The advantage of working
with bare ```Expression```s is that they are not bound to a specific model part or entity type, which means that they are completely interchangable. The following code snippet shows how to use bare ```Expression```s.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
model_part.CreateNewNode(3, 1.0, 1.0, 0.0)

prop = model_part.CreateNewProperties(1)
cond1 = model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
cond1.SetValue(Kratos.VELOCITY, Kratos.Array3([1,2,3]))

# now create the expression by reading non historical velocity:
exp = Kratos.Expression.IntegrationPointExpressionIO.Input(model_part, Kratos.VELOCITY, Kratos.Globals.DataLocation.Element).Execute()

# do some arithmetic operations
exp *= 2.0

# now write the expression value to model part as ACCELERATION in the historical container
Kratos.Expression.IntegrationPointExpressionIO.Output(model_part, Kratos.ACCELERATION, Kratos.Globals.DataLocation.Element).Execute(exp)
```
