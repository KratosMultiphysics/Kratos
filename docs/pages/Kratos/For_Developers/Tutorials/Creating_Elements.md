---
title: Creating the Element
keywords: 
tags: [Tutorial-Creating-the-Element.md]
sidebar: kratos_for_developers
summary: 
---

## Overview

As we wrote in the previous section, the files describing our element will be placed in the folder `custom_elements` , with names: `my_laplacian_element.h` and `my_laplacian_element.cpp`. The names of the subroutines are self-explanatory so their use is easy to understand. 

Note that since we have no source terms involving that involve the area of the elements (we have only puntual loads), no RHS term is calculated in the the subroutines. The Rigidity matrix is calculated in CalculateLocalSystem.

When the builder and solver is called in your strategy, it will loop through all the nodes and condition asking for their contribution to the system. This means that we need the following information:

* First of all, we need to know the contribution to the global system. In our case the 3x3 LHS matrix and the 3x1 RHS vector. Given the strategy we chose, to do we'll program this contribution in the CalculateLocalSystem. In this file we will perform all the necessary computations. We will return the rLeftHandSideMatrix (our LHS matrix) and the rRightHandSideVector (our RHS vector) . The rCurrentProcessInfo is an input used to gather general information, for example the `delta_time` in transient problems. In this particular problem we are not using it due to its simplicity.

* Besides returning the 3x3 matrix and the RHS vector, we must tell the builder to which degrees of freedom each of the lines of this matrix must add its contribution. In our case it is pretty simple, the first line corresponds to the `TEMPERATURE` DoF of the first node, the second node to the `TEMPERATURE` DoF of the second node and so on. This information is given to the builder and solver using both EquationIdVector and GetDofList.

In our case we will not use the CalculateRightHandSide, but the idea is basically the same as the CalculateLocalSystem.

## Header and source files

The application generator will provide us both header and source templates for the element. The header template contains the standard operations declaration, so we won't need to edit the header. We only need to edit the definition of mandatory operations in the source file:

* [EquationIdVector](Tutorial:-Creating-the-Element#equationidvector)
* [GetDofList](Tutorial:-Creating-the-Element#getdoflist)
* [CalculateLocalSystem](Tutorial:-Creating-the-Element#calculatelocalsystem)
* [Check](Tutorial:-Creating-the-Element#check)

The first lines inlcudes the license and the main authors. Do not forget to include the `geometry_utilities.h` and the `checks.h` files at the beginning of `my_laplacian_element.cpp`:

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
//  Main authors:    HERE_YOUR_NAME
//

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/variables.h"
#include "utilities/geometry_utilities.h"
#include "custom_elements/my_laplacian_element.h"
```


### EquationIdVector

```cpp
/**
 * this determines the elemental equation ID vector for all elemental DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void MyLaplacianElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const {
  unsigned int number_of_nodes = GetGeometry().PointsNumber();
  if(rResult.size() != number_of_nodes)
    rResult.resize(number_of_nodes,false);	

  for (unsigned int i=0;i<number_of_nodes;i++)
      rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
}
```


### GetDofList

```cpp
/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void MyLaplacianElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const {
  unsigned int number_of_nodes = GetGeometry().PointsNumber();
  if(rElementalDofList.size() != number_of_nodes)
    rElementalDofList.resize(number_of_nodes);

  for (unsigned int i=0;i<number_of_nodes;i++)
    rElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);
}
```


### CalculateLocalSystem

```cpp
/**
 * this is called during the assembling process in order
 * to calculate all elemental contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void MyLaplacianElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{

    KRATOS_TRY

    BoundedMatrix<double,3,2> DN_DX;  // Gradients matrix
    BoundedMatrix<double,2,2> D;      // Conductivity matrix
    D = ZeroMatrix(2,2); //initializing the matrix as zero
    array_1d<double,3> N; //dimension = number of nodes . Position of the gauss point
    array_1d<double,3> temp; //dimension = number of nodes . . since we are using a residualbased approach

    const unsigned int number_of_points = GetGeometry().size();

    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false); //resizing the system in case it does not have the right size

    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);

    // Getting data for the given geometry
    double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area); //asking for gradients and other info

    // Reading properties and conditions
    const double integrated_permittivity = area * GetProperties()[CONDUCTIVITY];
    D(0,0)=integrated_permittivity;
    D(1,1)=integrated_permittivity;

    // Main loop (one Gauss point)
    //const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    noalias(rLeftHandSideMatrix) = prod(DN_DX, Matrix(prod(D, trans(DN_DX))));  // Bt D B

    // Subtracting the dirichlet term
    // RHS -= LHS*DUMMY_UNKNOWNs
    for(unsigned int iii = 0; iii<number_of_points; iii++)
        temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,temp);

    KRATOS_CATCH("");
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
int MyLaplacianElement::Check(const ProcessInfo& rCurrentProcessInfo) const {

  KRATOS_TRY

      // Base class checks for positive Jacobian and Id > 0
      int ierr = Element::Check(rCurrentProcessInfo);
      if(ierr != 0) return ierr;

      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(TEMPERATURE)

      unsigned const int number_of_points = GetGeometry().size();  //added cornejo
      // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
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
