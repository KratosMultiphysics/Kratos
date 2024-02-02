---
title: Notes
keywords: 
tags: [notes, warnings, expression]
sidebar: kratos_expressions
summary: 
---

## Introduction

The use of ```Expression```, ```NodalExpression```, ```ConditionExpression``` or ```ElementExpression``` should be done carefully for the following reasons.

1. If it is possible to implement the arithmetics or opertations you do in ```C++``` level, then it is encouraged.
2. Using arithmetics operators should be done carefully. In the case if, an expression being updated by another vectorial expression using operators such as ```+=```, ```-=```, ```*=```, ```/=```, ```**=```, then it will keep both left hand side and right hand side double vectors in memory. If this is being done for lot of iterations, it can exceed available RAM and crash the execution. Therefore, it is a good practice to use ```Kratos.Expression.Utils.Collpase``` to collapse the lazy expression hierarchy.
3. Use of expressions has an overhead cost due to the use of runtime-polymorphism. So, if it is bottleneck, then transfer the expression to the ```C++```.
