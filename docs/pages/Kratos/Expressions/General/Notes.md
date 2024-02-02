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
2. Use arithmetic operators (`+=`, `-=`, `*=`, `/=`, `**=`, `+`, `-`, `*`, `/`) carefully. Operators keep their operands in memory, which means that if you chain operations for long enough, at some point you will run out of RAM. In such cases, `Kratos.Expression.Utils.Collpase` to forcefully evaluate an expression and discard the hierarchy of operands it stored. An example use case would be updating the same expression with new data in every time step of an analysis.
3. Use of expressions has an overhead cost due to the use of runtime-polymorphism. So, if it is bottleneck, then transfer the expression to the ```C++```.
