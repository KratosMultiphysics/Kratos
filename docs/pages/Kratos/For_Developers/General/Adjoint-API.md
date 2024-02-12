---
title: Adjoint API
keywords: 
tags: [Adjoint API]
sidebar: kratos_for_developers
summary: 
---

## Overview

This page describes the API for adjoint-based sensitivities in Kratos. It is the reference for developers creating new adjoint elements and response functions.

The main components for the adjoint and sensitivity calculations are the elements and conditions, the scheme and the response function. The elements are responsible for providing adjoint and sensitivity matrices which are gradients of the element residual vector. Everything related to the response function is provided by a separate class derived from the base class ResponseFunction. The response function is independent of the adjoint problem and can be used independently of strategies, schemes or builder and solvers. It uses the values stored in the adjoint variables to calculate sensitivities for a set of design variables given by the user. The scheme defines the left- and right-hand sides of the adjoint problem (transient or steady) based on the contributions from the adjoint element and response function and is independent of the sensitivity calculation. It is used in combination with the ResidualBasedLinearStrategy.

## Adjoint Variables

An adjoint variable is defined for each primal variable by prefixing the variable name with ADJOINT. For example, the adjoint variable of `DISPLACEMENT` is `ADJOINT_DISPLACEMENT`.

## Adjoint Elements and Conditions

Adjoint elements are created in the same way as regular elements. The adjoint element should be a separate class from the primal element. For applications with a relatively large number of linear elements like structural mechanics, the primal element type can be given as a template argument to a single generic adjoint element which stores the element internally and uses the element's public member functions to produce the adjoint matrices. Inheritance can also be used to minimize code duplication between primal and adjoint elements. The element functions and their meaning for the adjoint problem are given below.

### Solving the adjoint problem

The following functions are used in the solution of the adjoint problem:

```cpp
GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
```
Gets the list of adjoint dofs. This is the analog of the list of primal dofs.

```cpp
CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
```
Calculates adjoint matrix as the gradient of the negative residual with respect to the primal variable (e.g., DISPLACEMENT). For structural mechanics, this is the transpose of the stiffness matrix. This is zero if the primal variable is a rate (e.g., VELOCITY for fluid mechanics).

```cpp
CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
```
Calculates the adjoint matrix as the gradient of the negative residual with respect to first derivatives of the primal variable (e.g., VELOCITY). For structural mechanics, this is the transpose of the damping matrix. 

```cpp
CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
```
Calculates the adjoint matrix as the gradient of the negative residual with respect to second derivatives of the primal variable (e.g., ACCELERATION). This is normally the transpose of the mass matrix.

```cpp
GetValuesVector(Vector& values, int Step = 0)
```
Gets the vector of adjoint values for the primal variable (e.g., ADJOINT_DISPLACEMENT).

```cpp
GetFirstDerivativesVector(Vector& values, int Step = 0)
```
Gets the vector of adjoint values for first derivatives of primal variable (e.g., ADJOINT_VELOCITY).

```cpp
GetSecondDerivativesVector(Vector& values, int Step = 0)
```
Gets the vector of adjoint values for second derivatives of primal variable (e.g., ADJOINT_ACCELERATION).

### Calculating sensitivities

The following functions are used after the solution of the adjoint problem to calculate sensitivities:

```cpp
CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
```

Calculates the transpose of the gradient of the element residual with respect to the design variable. (These functions must be added to the element base class).

The sensitivities depend on the design variable which may be a `double` (e.g., THICKNESS) or `array_1d<double,3>` (e.g., NODAL_COORDINATES). The row dimension of the output is determined from the dimension of design variable (scalar or array_1d) and its storage type (nodal or elemental).

## Response Functions

The response function defines a scalar quantity of interest (e.g., drag), its gradients with respect to primal variables and their derivatives and its gradient with respect to design variables. The set of design variables is stored in the response function.

The following functions should be defined for each response function:

```cpp
CalculateValue()
```

Calculates the scalar value of the response function.

```cpp
CalculateGradient(const Element& rElem, const Matrix& rLHS, Vector& rOutput, ProcessInfo& rProcessInfo)
```

Calculates the gradient of the response function with respect to the primal variable (e.g., DISPLACEMENT). The rLHS matrix is the left hand side matrix of the corresponding adjoint equation.

```cpp
CalculateFirstDerivativesGradient(const Element& rElem, const Matrix& rLHS, Vector& rOutput, ProcessInfo& rProcessInfo)
```

Calculates the gradient of the response function with respect to first derivatives of the primal variable (e.g., VELOCITY).

```cpp
CalculateSecondDerivativesGradient(const Element& rElem, const Matrix& rLHS, Vector& rOutput, ProcessInfo& rProcessInfo)
```

Calculates the gradient of the response function with respect to second derivatives of the primal variable (e.g., ACCELERATION).

```cpp
UpdateSensitivities()
```

Updates the sensitivities for the design variables based on the most recent adjoint solution. The design variables are specified in the json parameters of the response function.