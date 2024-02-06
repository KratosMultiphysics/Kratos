---
title: Domain Size Expression IO
keywords: 
tags: [domain size, condition, element, expression io]
sidebar: kratos_expressions
summary: 
---

## Introduction

```DomainSizeExpressionIO``` computes the domain sizes (length, surface or volume, depending of the entity's dimension) of each entitiy (such as ```Condition``` or ```Element```) in the provided `ModelPart`. This only supports ```ConditionExpression``` or ```ElementExression```. Domain sizes are computed by accessing underlying geometry's ```Geometry::DomainSize``` method.

1. If the ```Condition``` or ```Element``` as a geometry representing a line, then this will compute an expression with line length for each entitiy.
2. If the ```Condition``` or ```Element``` as a geometry representing a surface, then this will compute an expression with surface area for each entitiy.
3. If the ```Condition``` or ```Element``` as a geometry representing a volume, then this will compute an expression with volume for each entitiy.

The resulting expression is always a scalar expression.

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
(1,)
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
(DoubleVec[]*2)
```
