---
title: Expression Types
keywords: 
tags: [Expression_Types.md]
sidebar: kratos_expressions
summary: 
---

## Introduction

There are four types of expressions.
* `Expression`
* `NodalExpression`
* `ConditionExpression`
* `ElementExpression`

## `Expression`

`Expression` is the lowest level class that everything else is built on top of. It models an array of numeric values (scalars, vectors, matrices, etc.), supports arithmetic operations, but is not related to any container in Kratos. To create an `Expression` or write an existing one, you can use one of the derived classes of `ExpressionIO`.

TODO: add an example (@matekelemen) - requires exposing the `Input` and `Output` nested classes of `VariableExpressionIO` to python.


## Nodal Expression

Nodal expressions are used to store data related to nodal quantities. They may be historical or non-historical. Following code snippet illustrtes creating a nodal expression
and reading ```VELOCITY``` from non-historical nodal container.
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
```

## Condition Expression

Condition expressions are used to store data related to condition quantities. Following code snippet illustrtes creating a condition expression
and reading ```VELOCITY``` from the conditions' container.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)

prop = model_part.CreateNewProperties(1)
cond_1 = model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], prop)
cond_1.SetValue(Kratos.VELOCITY, Kratos.Array3([1,2,3]))

# now create the expression:
condition_expression = Kratos.Expression.ConditionExpression(model_part)

# now read the VELOCITY from the condition data value container
Kratos.Expression.VariableExpressionIO.Read(condition_expression, Kratos.VELOCITY)
```

## Element Expression

Element expressions are used to store data related to element quantities. Following code snippet illustrtes creating a element expression
and reading ```VELOCITY``` from the elements' container.
```python
import KratosMultiphysics as Kratos
model = Kratos.Model()
model_part = model.CreateModelPart("test")
node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node_2 = model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
node_3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)

prop = model_part.CreateNewProperties(1)
elem_1 = model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
elem_1.SetValue(Kratos.VELOCITY, Kratos.Array3([1,2,3]))

# now create the expression:
element_expression = Kratos.Expression.ElementExpression(model_part)

# now read the VELOCITY from the condition data value container
Kratos.Expression.VariableExpressionIO.Read(element_expression, Kratos.VELOCITY)
```
