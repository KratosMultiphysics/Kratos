---
title: Overview
keywords: 
tags: [overview.md]
sidebar: kratos_expressions
summary: 
---
## Introduction

An expression represents a mechanism where a given mathematical expression is evaluated when actual value of that expression in required. One can compare these
expressions to lazy expressions (expression templates) in [Eigen](https://eigen.tuxfamily.org/dox/TopicLazyEvaluation.html) or in [Boost](https://www.boost.org/doc/libs/1_82_0/doc/html/boost_yap/manual.html). The main difference being, Kratos expressions are **dynamic** lazy expressions (evaluated on demand, not automatically) whose type is evaluated at runtime, whereas Eigen and Boost has **static** lazy expressions whose type is known at compile-time. Few of the main advantages of these lazy expressions are:
* May reduce memory footprint.
* May reduce computational cost.

KratosMultiphysics stores its data in different containers depending on what kind of geometry that data is associated with (node, element, condition, etc.) and what that data represents. Each container requires a ```Kratos::Variable``` to access its underlying data. This requires defining the variable beforehand, and then using it to store and retrieve its value. Expressions allow users/developers to access and manipulate the data within these containers regardless of what `Kratos::Variable` they belong to. All expressions support **shared memory** and **distributed memory** parallelism. Though implemented in C++, expressions are meant to be used mainly in python.

## Supported data containers

The following data containers are supported:
* NodalExpression - data related to nodes (with historical or non-historical variables).
* ConditionExpression - data related to conditions.
* ElementExpression - data related to elements.

## Supported data types

Expressions are compatible with most main variable types in Kratos:

* ```int```
* ```double```
* ```array_1d<double, 3>```
* ```array_1d<double, 4>```
* ```array_1d<double, 6>```
* ```array_1d<double, 9>```
* ```Vector```
* ```Matrix```

## Use cases

Expressions are useful in following scenarios:
1. Transferring data between Kratos and other libraries - third-party libraries may be developed in other languages like Fortran as well.
2. Transferring data between Kratos and numpy/scipy.
3. Storing and fetching Kratos data in HDF5.
4. Non performance-critical intermediate calculations in python.
5. Transerring data back and forth within Kratos without variables such as in OptimizationApplication.
6. Quick and convenient prototyping of research code in python while using **shared memory** and **distributed memory** paralellism.





