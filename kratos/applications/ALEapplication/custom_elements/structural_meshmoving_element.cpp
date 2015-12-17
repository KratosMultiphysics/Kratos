/*
 ==============================================================================
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 Version 1.0 (Released on march 05, 2007).

 Copyright 2007
 Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
 pooyan@cimne.upc.edu
 rrossi@cimne.upc.edu
 janosch.stascheit@rub.de
 nagel@sd.rub.de
 - CIMNE (International Center for Numerical Methods in Engineering),
 Gran Capita' s/n, 08034 Barcelona, Spain
 - Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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

 ==============================================================================
 */

/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: Andreas.Mini@tum.de $
 *  Date:                $Date: November 2015 $
 *  Revision:            $Revision: 1.4 $
 * ***************************************************************************/

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/structural_meshmoving_element.h"
#include "ale_application.h"

namespace Kratos {

StructuralMeshMovingElement::StructuralMeshMovingElement(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {
  mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

}

StructuralMeshMovingElement::StructuralMeshMovingElement(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {
  mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

void StructuralMeshMovingElement::GetDisplacementValues(VectorType& rValues,
                                                        const int Step) {
  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType LocalSize = NumNodes * dimension;

  if (rValues.size() != LocalSize)
    rValues.resize(LocalSize, false);

  if (dimension == 2) {

    SizeType Index = 0;

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {

      rValues[Index++] = rGeom[iNode].FastGetSolutionStepValue(DISPLACEMENT_X,
                                                               Step);
      rValues[Index++] = rGeom[iNode].FastGetSolutionStepValue(DISPLACEMENT_Y,
                                                               Step);
    }

  } else if (dimension == 3) {

    SizeType Index = 0;

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {

      rValues[Index++] = rGeom[iNode].FastGetSolutionStepValue(DISPLACEMENT_X,
                                                               Step);
      rValues[Index++] = rGeom[iNode].FastGetSolutionStepValue(DISPLACEMENT_Y,
                                                               Step);
      rValues[Index++] = rGeom[iNode].FastGetSolutionStepValue(DISPLACEMENT_Z,
                                                               Step);
    }
  }
}

void StructuralMeshMovingElement::Initialize() {

  const GeometryType::IntegrationPointsArrayType& integration_points =
      GetGeometry().IntegrationPoints(mThisIntegrationMethod);

  mInvJ0.resize(integration_points.size());
  mDetJ0.resize(integration_points.size(), false);

  this->GetGeometry().Jacobian(mJ0, mThisIntegrationMethod);
  this->GetGeometry().InverseOfJacobian(mInvJ0, mThisIntegrationMethod);
  this->GetGeometry().DeterminantOfJacobian(mDetJ0);
}

StructuralMeshMovingElement::MatrixType StructuralMeshMovingElement::SetAndModifyConstitutiveLaw(
    const int &dimension) {

  KRATOS_TRY;

  double YoungsModulus = 200000;
  const double PoissonCoefficient = 0.3;

  // Stiffening of elements using Jacobian determinants and exponent between 0.0 and 2.0
  const double J0 = 100;  // Factor influences how far the displacement spreads into the fluid mesh
  const double Xi = 1.5;  // 1.5 Exponent influences stiffening of smaller elements; 0 = no stiffening
  const double DetJMag = std::abs(mDetJ0[0]);
  const double Quotient = J0 / DetJMag;
  YoungsModulus *= (DetJMag * pow(Quotient, Xi));

  double Lambda = (YoungsModulus * PoissonCoefficient)
      / ((1 + PoissonCoefficient) * (1 - 2 * PoissonCoefficient));
  double Mue = YoungsModulus / (2 * (1 - PoissonCoefficient));

  MatrixType ConstitutiveMatrix;

  if (dimension == 2) {
    ConstitutiveMatrix = ZeroMatrix(3, 3);

    ConstitutiveMatrix(0, 0) = Lambda + 2 * Mue;
    ConstitutiveMatrix(1, 1) = ConstitutiveMatrix(0, 0);
    ConstitutiveMatrix(2, 2) = Mue;
    ConstitutiveMatrix(0, 1) = Lambda;
    ConstitutiveMatrix(1, 0) = Lambda;
  }

  else if (dimension == 3) {

    ConstitutiveMatrix = ZeroMatrix(6, 6);

    double ConstitutiveMatrixComponent01 = (YoungsModulus
        * (1 - PoissonCoefficient)
        / ((1 + PoissonCoefficient) * (1 - 2 * PoissonCoefficient)));

    double ConstitutiveMatrixComponent02 = YoungsModulus
        * (1 - PoissonCoefficient) / (1 + PoissonCoefficient);

    double ConstitutiveMatrixComponent03 = YoungsModulus * PoissonCoefficient
        / ((1 + PoissonCoefficient) * (1 - 2 * PoissonCoefficient));

    ConstitutiveMatrix(0, 0) = ConstitutiveMatrixComponent01;
    ConstitutiveMatrix(1, 1) = ConstitutiveMatrixComponent01;
    ConstitutiveMatrix(2, 2) = ConstitutiveMatrixComponent01;
    ConstitutiveMatrix(3, 3) = ConstitutiveMatrixComponent02;
    ConstitutiveMatrix(4, 4) = ConstitutiveMatrixComponent02;
    ConstitutiveMatrix(5, 5) = ConstitutiveMatrixComponent02;

    ConstitutiveMatrix(0, 1) = ConstitutiveMatrixComponent03;
    ConstitutiveMatrix(1, 0) = ConstitutiveMatrixComponent03;
    ConstitutiveMatrix(0, 2) = ConstitutiveMatrixComponent03;
    ConstitutiveMatrix(2, 0) = ConstitutiveMatrixComponent03;
    ConstitutiveMatrix(1, 2) = ConstitutiveMatrixComponent03;
    ConstitutiveMatrix(2, 1) = ConstitutiveMatrixComponent03;

  }

  return ConstitutiveMatrix;

  KRATOS_CATCH("");
}

StructuralMeshMovingElement::MatrixType StructuralMeshMovingElement::CalculateBMatrix(
    const int &dimension) {

  KRATOS_TRY;

  GeometryType::JacobiansType j;

  GetGeometry().Jacobian(j, mThisIntegrationMethod);

  Matrix F = prod(j[0], mInvJ0[0]);

  GeometryType::ShapeFunctionsGradientsType DN_De = this->GetGeometry()
      .ShapeFunctionsLocalGradients(mThisIntegrationMethod);

  Matrix DN_DX = prod(DN_De[0], mInvJ0[0]);

  const SizeType NumNodes = this->GetGeometry().PointsNumber();

  MatrixType B;

  if (dimension == 2) {

    B = ZeroMatrix(3, NumNodes * 2);

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {

      SizeType index = 2 * iNode;

      B(0, index + 0) = F(0, 0) * DN_DX(iNode, 0);
      B(0, index + 1) = F(1, 0) * DN_DX(iNode, 0);
      B(1, index + 0) = F(0, 1) * DN_DX(iNode, 1);
      B(1, index + 1) = F(1, 1) * DN_DX(iNode, 1);
      B(2, index + 0) = F(0, 0) * DN_DX(iNode, 1) + F(0, 1) * DN_DX(iNode, 0);
      B(2, index + 1) = F(1, 0) * DN_DX(iNode, 1) + F(1, 1) * DN_DX(iNode, 0);
    }
  }

  else if (dimension == 3) {

    B = ZeroMatrix(6, NumNodes * 3);

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {

      SizeType index = 3 * iNode;

      B(0, index + 0) = F(0, 0) * DN_DX(iNode, 0);
      B(0, index + 1) = F(1, 0) * DN_DX(iNode, 0);
      B(0, index + 2) = F(2, 0) * DN_DX(iNode, 0);
      B(1, index + 0) = F(0, 1) * DN_DX(iNode, 1);
      B(1, index + 1) = F(1, 1) * DN_DX(iNode, 1);
      B(1, index + 2) = F(2, 1) * DN_DX(iNode, 1);
      B(2, index + 0) = F(0, 2) * DN_DX(iNode, 2);
      B(2, index + 1) = F(1, 2) * DN_DX(iNode, 2);
      B(2, index + 2) = F(2, 2) * DN_DX(iNode, 2);
      B(3, index + 0) = F(0, 0) * DN_DX(iNode, 1) + F(0, 1) * DN_DX(iNode, 0);
      B(3, index + 1) = F(1, 0) * DN_DX(iNode, 1) + F(1, 1) * DN_DX(iNode, 0);
      B(3, index + 2) = F(2, 0) * DN_DX(iNode, 1) + F(2, 1) * DN_DX(iNode, 0);
      B(4, index + 0) = F(0, 1) * DN_DX(iNode, 2) + F(0, 2) * DN_DX(iNode, 1);
      B(4, index + 1) = F(1, 1) * DN_DX(iNode, 2) + F(1, 2) * DN_DX(iNode, 1);
      B(4, index + 2) = F(2, 1) * DN_DX(iNode, 2) + F(2, 2) * DN_DX(iNode, 1);
      B(5, index + 0) = F(0, 2) * DN_DX(iNode, 0) + F(0, 0) * DN_DX(iNode, 2);
      B(5, index + 1) = F(1, 2) * DN_DX(iNode, 0) + F(1, 0) * DN_DX(iNode, 2);
      B(5, index + 2) = F(2, 2) * DN_DX(iNode, 0) + F(2, 0) * DN_DX(iNode, 2);
    }
  }

  return B;

  KRATOS_CATCH("");
}

void StructuralMeshMovingElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {

  KRATOS_TRY;

  const SizeType NumNodes = this->GetGeometry().PointsNumber();
  const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
  const SizeType LocalSize = NumNodes * dimension;

  MatrixType ConstitutiveMatrix = SetAndModifyConstitutiveLaw(dimension);

  MatrixType B = CalculateBMatrix(dimension);

  // Compute LHS
  if (rLeftHandSideMatrix.size1() != LocalSize)
    rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

  MatrixType intermediateMatrix = prod(trans(B), ConstitutiveMatrix);
  noalias(rLeftHandSideMatrix) = prod(intermediateMatrix, B);

  // Compute RHS
  if (rRightHandSideVector.size() != LocalSize)
    rRightHandSideVector.resize(LocalSize, false);

  VectorType LastValues = ZeroVector(LocalSize);
  this->GetDisplacementValues(LastValues, 0);
  noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, LastValues);

  KRATOS_CATCH("");
}

void StructuralMeshMovingElement::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {

  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.size();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType LocalSize = NumNodes * dimension;

  if (rResult.size() != LocalSize)
    rResult.resize(LocalSize, false);

  if (dimension == 2) {

    SizeType Index = 0;

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {

      rResult[Index++] = rGeom[iNode].GetDof(DISPLACEMENT_X).EquationId();
      rResult[Index++] = rGeom[iNode].GetDof(DISPLACEMENT_Y).EquationId();
    }
  } else if (dimension == 3) {

    SizeType Index = 0;

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {

      rResult[Index++] = rGeom[iNode].GetDof(DISPLACEMENT_X).EquationId();
      rResult[Index++] = rGeom[iNode].GetDof(DISPLACEMENT_Y).EquationId();
      rResult[Index++] = rGeom[iNode].GetDof(DISPLACEMENT_Z).EquationId();
    }
  }
}

void StructuralMeshMovingElement::GetDofList(DofsVectorType& rElementalDofList,
                                             ProcessInfo& rCurrentProcessInfo) {

  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.size();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType LocalSize = NumNodes * dimension;

  if (rElementalDofList.size() != LocalSize)
    rElementalDofList.resize(LocalSize);

  if (dimension == 2) {

    SizeType LocalIndex = 0;

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {

      rElementalDofList[LocalIndex++] = rGeom[iNode].pGetDof(DISPLACEMENT_X);
      rElementalDofList[LocalIndex++] = rGeom[iNode].pGetDof(DISPLACEMENT_Y);
    }
  } else if (dimension == 3) {

    SizeType LocalIndex = 0;

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode) {
      rElementalDofList[LocalIndex++] = rGeom[iNode].pGetDof(DISPLACEMENT_X);
      rElementalDofList[LocalIndex++] = rGeom[iNode].pGetDof(DISPLACEMENT_Y);
      rElementalDofList[LocalIndex++] = rGeom[iNode].pGetDof(DISPLACEMENT_Z);
    }
  }
}

}  // Namespace Kratos
