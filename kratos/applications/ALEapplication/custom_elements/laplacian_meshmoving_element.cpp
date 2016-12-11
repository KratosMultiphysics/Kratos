// ==============================================================================
/*
 KratosALEApllication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Pooyan Dadvand, Riccardo Rossi, Andreas Winterstein
                     pooyan@cimne.upc.edu
                     rrossi@cimne.upc.edu
                     a.winterstein@tum.de
- CIMNE (International Center for Numerical Methods in Engineering),
  Gran Capita' s/n, 08034 Barcelona, Spain
- Chair of Structural Analysis, Technical University of Munich
  Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================

/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: A.Winterstein@tum.de $
 *  Date:                $Date: June 2016 $
 *  Revision:            $Revision: 1.5 $
 * ***************************************************************************/

// System includes
// External includes
// Project includes
#include "includes/define.h"
#include "ale_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "custom_elements/laplacian_meshmoving_element.h"

namespace Kratos {

LaplacianMeshMovingElement::LaplacianMeshMovingElement(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {
  mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

LaplacianMeshMovingElement::LaplacianMeshMovingElement(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {
  mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

Element::Pointer LaplacianMeshMovingElement::Create(
    IndexType NewId, NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const {

  return Element::Pointer(
      new LaplacianMeshMovingElement(NewId, GetGeometry().Create(rThisNodes),
                                     pProperties));
}

void LaplacianMeshMovingElement::Initialize() {

  KRATOS_TRY;

  const GeometryType::IntegrationPointsArrayType& integration_points =
      GetGeometry().IntegrationPoints(mThisIntegrationMethod);

  mJ0 = GetGeometry().Jacobian(mJ0, mThisIntegrationMethod);
  mInvJ0.resize(integration_points.size());
  mDetJ0.resize(integration_points.size(), false);
  mTotalDomainInitialSize = 0.00;

  for (unsigned int PointNumber = 0; PointNumber < integration_points.size();
      ++PointNumber) {
    //getting informations for integration
    double IntegrationWeight = integration_points[PointNumber].Weight();

    //calculating and storing inverse of the jacobian and the parameters needed
    MathUtils<double>::InvertMatrix(mJ0[PointNumber], mInvJ0[PointNumber],
                                    mDetJ0[PointNumber]);

    //calculating the total area
    mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
  }

  KRATOS_CATCH("");

}

LaplacianMeshMovingElement::MatrixType LaplacianMeshMovingElement::CalculateDerivatives(
    const int& rdimension, const double& rPointNumber) {

  KRATOS_TRY;

  GeometryType::ShapeFunctionsGradientsType DN_De = this->GetGeometry()
      .ShapeFunctionsLocalGradients(mThisIntegrationMethod);

  Matrix DN_DX = prod(DN_De[rPointNumber], mInvJ0[rPointNumber]);

  return DN_DX;

  KRATOS_CATCH("");

}

void LaplacianMeshMovingElement::CalculateDeltaPosition(
    VectorType& IntermediateDisplacements, ProcessInfo& rCurrentProcessInfo) {

  KRATOS_TRY;

  unsigned int ComponentIndex = rCurrentProcessInfo[FRACTIONAL_STEP] - 1;

  const SizeType NumNodes = this->GetGeometry().PointsNumber();

  for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {

    const VectorType& displacement = GetGeometry()[iNode]
        .FastGetSolutionStepValue(MESH_DISPLACEMENT, 0)
        - GetGeometry()[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
    IntermediateDisplacements[iNode] = displacement[ComponentIndex];
  }

  KRATOS_CATCH("");
}

void LaplacianMeshMovingElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {

  KRATOS_TRY;

  const SizeType NumNodes = this->GetGeometry().PointsNumber();
  const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

  const GeometryType::IntegrationPointsArrayType& integration_points =
      GetGeometry().IntegrationPoints(mThisIntegrationMethod);

  if (rLeftHandSideMatrix.size1() != NumNodes)
    rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

  if (rRightHandSideVector.size() != NumNodes)
    rRightHandSideVector.resize(NumNodes, false);

  noalias(rLeftHandSideMatrix) = ZeroMatrix(NumNodes, NumNodes);
  noalias(rRightHandSideVector) = ZeroVector(NumNodes);

  for (unsigned int PointNumber = 0; PointNumber < integration_points.size();
      ++PointNumber) {

    double IntegrationWeight = integration_points[PointNumber].Weight();

    //Compute LHS
    MatrixType DN_DX = CalculateDerivatives(dimension, PointNumber);
    noalias(rLeftHandSideMatrix) += prod((IntegrationWeight * DN_DX),
                                         trans(DN_DX));

    VectorType IntermediateDisplacements(NumNodes);

    CalculateDeltaPosition(IntermediateDisplacements, rCurrentProcessInfo);

    //Compute RHS
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, IntermediateDisplacements);

  }KRATOS_CATCH("");

}

void LaplacianMeshMovingElement::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {

  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.size();
  unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType LocalSize = NumNodes * dimension;

  if (rResult.size() != LocalSize)
    rResult.resize(LocalSize, false);

  for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {

    SizeType Index = iNode * dimension;
    rResult[Index] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_X).EquationId();
    rResult[Index + 1] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_Y).EquationId();
    if (dimension == 3)
      rResult[Index + 2] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_Z).EquationId();
  }
}

void LaplacianMeshMovingElement::GetDofList(DofsVectorType& rElementalDofList,
                                            ProcessInfo& rCurrentProcessInfo) {

  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.size();
  unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType LocalSize = NumNodes * dimension;

  if (rElementalDofList.size() != LocalSize)
    rElementalDofList.resize(LocalSize);

  for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {

    SizeType Index = iNode * dimension;

    rElementalDofList[Index] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_X);
    rElementalDofList[Index + 1] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_Y);
    if (dimension == 3)
      rElementalDofList[Index + 2] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_Z);

  }
}

}  // Namespace Kratos
