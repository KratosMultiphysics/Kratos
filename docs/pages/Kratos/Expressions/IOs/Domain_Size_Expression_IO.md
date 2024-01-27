---
title: Domain Size Expression IO
keywords:
tags: [domain size, condition, element, expression io]
sidebar: kratos_expressions
summary:
---

## Introduction

```DomainSizeExpressionIO``` computes the domain sizes of each entitiy (such as ```Condition``` or ```Element```). This only supports ```ConditionExpression``` or ```ElementExression```. The domain sizes are computed by accessing underlying geometry's ```Geometry::DomainSize``` method of the entity.

1. If the ```Condition``` or ```Element``` as a geometry representing a line, then this will compute an expression with line length for each entitiy.
2. If the ```Condition``` or ```Element``` as a geometry representing a surface, then this will compute an expression with surface area for each entitiy.
3. If the ```Condition``` or ```Element``` as a geometry representing a volume, then this will compute an expression with volume for each entitiy.

The resulting expression will be always a scalar expression.

## Reading domain sizes
Following code snippet shows how to read domain sizes to an expression
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
model_part.CreateNewNode(3, 1.0, 1.0, 0.0)

prop = model_part.CreateNewProperties(1)
model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

# now create the expression:
element_expression = Kratos.Expression.ElementExpression(model_part)

# now read the element sizes (in this case triangle surface area)
Kratos.Expression.DomainSizeExpressionIO.Read(element_expression)

shape = element_expression.Evaluate().shape
print(shape)
```

## Using expressions without the model parts
The ```ConditionExpression``` and ```ElementExpression``` has an expression which can be directly used if required. The advantage of working
with the ```Expression``` directely is, then it is not bound to a model part of a ```DataValueContainer```. Hence, these expressions can be interchanged if required in
advanced use cases. Following code snippet shows how to use bare ```Expressions```.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
model_part.CreateNewNode(3, 1.0, 1.0, 0.0)

prop = model_part.CreateNewProperties(1)
model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

# now create the expression by reading element domain sizes:
exp = Kratos.Expression.DomainSizeExpressionIO.Input(model_part, Kratos.Globals.DataLocation.Element).Execute()

# do some arithmetic operations
exp *= 2.0

print(exp)
```
