---
title: Creating the Conditions
keywords: 
tags: [Tutorial-Creating-the-Conditions.md]
sidebar: kratos_for_developers
summary: 
---

## Overview
 
What you need to understand of these lines

When the builder and solver is called in your strategy, it will loop through all the nodes and condition asking for their contribution to the system. This means that we need the following information:

* First of all, we need to know the contribution to the global system. In our case the 1x1 LHS matrix and the 1x1 RHS vector. Given the strategy we chose, to do we'll program this contribution in the **CalculateLocalSystem**. In this file we will perform all the necessary computations. We will return the rLeftHandSideMatrix (our LHS matrix) and the rRightHandSideVector (our RHS vector) . The rCurrentProcessInfo is an input used to gather general information, for example the `delta_time` in transcient problems. In this particular problem we are not using it due to its simplicity. Since it is only a load, we will not modify the LHS matrix but only the RHS vector.

* Besides returning the 1x1 matrix and the RHS vector, we must tell the builder to which degrees of freedom each of the lines of this matrix must add its contribution. In our case it is pretty simple, the first line corresponds to the `TEMPERATURE` DoF of the first node, the second node to the `TEMPERATURE` DoF of the second node and so on. This information is given to the builder and solver using both **EquationIdVector** and **GetDofList**.

In this particular case, since our condition is only linked to one Degree of Freedom, it has a 1x1 matrix and a 1x1 RHS vector. If we were using a line condition (2 nodes), our matrix would be 2x2 and the RHS 2x1.

In our case we will not use the CalculateRightHandSide, but the idea is basically the same as the CalculateLocalSystem.

##  Creating a Custom Condition to incorporate the thermal loads

The process of creating a custom condition is similar to creating an element, and the application generator provides us both header and source templates.

`POINT_HEAT_SOURCE` is our variable containing the external thermal loads to our problem. KRATOS will call the condition and assign it node by node to the system of equations. Since this is standard load, we only have to add the contribution to the right hand side (RHS) without any calculation. In pseudo code, our task is simply: 

```cpp
RHS (node) + = POINT_HEAT_SOURCE (node)
```

Still, note that the structure of a condition, as well as the one in elements, includes both the RHS and the left hand side matrix (LHS). This gives us flexibility and thus we can incorporate operations that also affect the matrix. Bear in mind that conditions are a class called by the **solver**, so we will edit the mandatory operations:

* [EquationIdVector](Tutorial:-Creating-the-Conditions#equationidvector)
* [GetDofList](Tutorial:-Creating-the-Conditions#getdoflist)
* [CalculateLocalSystem](Tutorial:-Creating-the-Conditions#calculatelocalsystem)
* [Check](Tutorial:-Creating-the-Conditions#check)

Just as in the elements, to incorporate this condition we have to create it and register it in Kratos. Once registered, the solver will be the one calling it and assembling the LHS and RHS, so there's no need to call it in the problemtype.

Here are the first lines of `point_source_condition.cpp`:

```cpp
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    , KratosAppGenerator
//

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/checks.h"
#include "includes/variables.h"
#include "custom_conditions/point_source_condition.h"
```


### EquationIdVector

```cpp
/**
 * this determines the condition equation ID vector for all conditional
 * DOFs
 * @param rResult: the condition equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes, false);      // False says not to preserve existing storage!!

    for (unsigned int i = 0; i < number_of_nodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
}
```	

### GetDofList

```cpp
/**
 * determines the condition equation list of DOFs
 * @param ConditionDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rConditionDofList.size() != number_of_nodes)
        rConditionDofList.resize(number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; i++)
        rConditionDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);
}
```

### CalculateLocalSystem

```cpp
/**
 * this is called during the assembling process in order
 * to calculate all condition contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rRightHandSideVector: the condition right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rLeftHandSideMatrix.size1() != 1)
        rLeftHandSideMatrix.resize(1,1,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(1,1);

    if(rRightHandSideVector.size() != 1)
        rRightHandSideVector.resize(1,false);
    double load = GetGeometry()[0].GetSolutionStepValue(POINT_HEAT_SOURCE);
    rRightHandSideVector[0] = load;
    KRATOS_CATCH("")
}
```

### Check

```cpp
/**
 * This method provides the place to perform checks on the completeness of the input
 * and the compatibility with the problem options as well as the contitutive laws selected
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
int PointSourceCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Condition::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(TEMPERATURE)
    KRATOS_CHECK_VARIABLE_KEY(POINT_HEAT_SOURCE)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const unsigned int number_of_points = GetGeometry().size();
    for ( unsigned int i = 0; i < number_of_points; i++ )
    {
        auto &rnode = this->GetGeometry()[i];
        KRATOS_EXPECT_VARIABLE_IN_NODAL_DATA(TEMPERATURE,rnode)
        KRATOS_EXPECT_DOF_IN_NODE(TEMPERATURE,rnode)
    }

    return ierr;

    KRATOS_CATCH("");
}
```
