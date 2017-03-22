// =============================================================================
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
#include "custom_elements/structural_meshmoving_element.h"
#include "ale_application.h"

namespace Kratos
{
StructuralMeshMovingElement::StructuralMeshMovingElement(IndexType NewId,
                                                         GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

StructuralMeshMovingElement::StructuralMeshMovingElement(IndexType NewId,
                                                         GeometryType::Pointer pGeometry,
                                                         PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

Element::Pointer StructuralMeshMovingElement::Create(IndexType NewId,
                                                     NodesArrayType const& rThisNodes,
                                                     PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = this->GetGeometry();

    return BaseType::Pointer(new StructuralMeshMovingElement(
        NewId, rGeom.Create(rThisNodes), pProperties));
}

Element::Pointer StructuralMeshMovingElement::Create(IndexType NewId,
                                                     GeometryType::Pointer pGeom,
                                                     PropertiesType::Pointer pProperties) const
{
    return BaseType::Pointer(new StructuralMeshMovingElement(NewId, pGeom, pProperties));
}

void StructuralMeshMovingElement::GetDisplacementValues(VectorType& rValues, const int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    if (dimension == 2)
    {
        SizeType Index = 0;

        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            rValues[Index++] =
                rGeom[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT_X, Step);
            rValues[Index++] =
                rGeom[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT_Y, Step);
        }
    }
    else if (dimension == 3)
    {
        SizeType Index = 0;

        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            rValues[Index++] =
                rGeom[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT_X, Step);
            rValues[Index++] =
                rGeom[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT_Y, Step);
            rValues[Index++] =
                rGeom[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT_Z, Step);
        }
    }
}

void StructuralMeshMovingElement::Initialize()
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);

    mInvJ0.resize(integration_points.size());
    mDetJ0.resize(integration_points.size(), false);
    mTotalDomainInitialSize = 0.00;

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    mLocalSize = NumNodes * dimension;

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        // getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();

        // calculating and storing inverse of the jacobian and the parameters
        // needed
        MathUtils<double>::InvertMatrix(
            J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber]);

        // calculating the total area
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }

    KRATOS_CATCH("");
}

StructuralMeshMovingElement::MatrixType StructuralMeshMovingElement::SetAndModifyConstitutiveLaw(
    const int& dimension, const double& rPointNumber)
{
    KRATOS_TRY;

    double YoungsModulus = 200000;
    const double PoissonCoefficient = 0.3;

    // double detJ = mDetJ0[rPointNumber] * IntegrationWeight;

    VectorType DetJ;
    this->GetGeometry().DeterminantOfJacobian(DetJ);

    // Stiffening of elements using Jacobian determinants and exponent between
    // 0.0 and 2.0
    const double J0 = 100; // Factor influences how far the displacement spreads
                           // into the fluid mesh
    const double Xi = 1.5; // 1.5 Exponent influences stiffening of smaller
                           // elements; 0 = no stiffening
    // double DetJMag = detJ;
    double DetJMag = mTotalDomainInitialSize; // std::abs(DetJ[0]);
    const double Quotient = J0 / DetJMag;
    YoungsModulus *= (DetJMag * pow(Quotient, Xi));

    double Lambda = (YoungsModulus * PoissonCoefficient) /
                    ((1 + PoissonCoefficient) * (1 - 2 * PoissonCoefficient));
    double Mue = YoungsModulus / (2 * (1 - PoissonCoefficient));

    MatrixType ConstitutiveMatrix;

    if (dimension == 2)
    {
        ConstitutiveMatrix = ZeroMatrix(3, 3);

        ConstitutiveMatrix(0, 0) = Lambda + 2 * Mue;
        ConstitutiveMatrix(1, 1) = ConstitutiveMatrix(0, 0);
        ConstitutiveMatrix(2, 2) = Mue;
        ConstitutiveMatrix(0, 1) = Lambda;
        ConstitutiveMatrix(1, 0) = Lambda;
    }

    else if (dimension == 3)
    {
        ConstitutiveMatrix = ZeroMatrix(6, 6);

        double ConstitutiveMatrixComponent01 =
            (YoungsModulus * (1 - PoissonCoefficient) /
             ((1 + PoissonCoefficient) * (1 - 2 * PoissonCoefficient)));

        double ConstitutiveMatrixComponent02 =
            YoungsModulus * (1 - PoissonCoefficient) / (1 + PoissonCoefficient);

        double ConstitutiveMatrixComponent03 =
            YoungsModulus * PoissonCoefficient /
            ((1 + PoissonCoefficient) * (1 - 2 * PoissonCoefficient));

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
    const int& dimension, const double& rPointNumber)
{
    KRATOS_TRY;

    GeometryType::JacobiansType j;

    GetGeometry().Jacobian(j, mThisIntegrationMethod);

    Matrix F = prod(j[rPointNumber], mInvJ0[rPointNumber]);

    GeometryType::ShapeFunctionsGradientsType DN_De =
        this->GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    Matrix DN_DX = prod(DN_De[rPointNumber], mInvJ0[rPointNumber]);

    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    MatrixType B;

    if (dimension == 2)
    {
        B = ZeroMatrix(3, NumNodes * 2);

        // SizeType Index = 0;

        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType Index = 2 * iNode;

            B(0, Index + 0) = DN_DX(iNode, 0);
            B(0, Index + 1) = 0.0;
            B(1, Index + 0) = 0.0;
            B(1, Index + 1) = DN_DX(iNode, 1);
            B(2, Index + 0) = DN_DX(iNode, 1);
            B(2, Index + 1) = DN_DX(iNode, 0);
        }
    }

    else if (dimension == 3)
    {
        B = ZeroMatrix(6, NumNodes * 3);

        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType index = 3 * iNode;

            B(0, index + 0) = DN_DX(iNode, 0);
            B(1, index + 1) = DN_DX(iNode, 1);
            B(2, index + 2) = DN_DX(iNode, 2);

            B(3, index + 0) = DN_DX(iNode, 1);
            B(3, index + 1) = DN_DX(iNode, 0);

            B(4, index + 1) = DN_DX(iNode, 2);
            B(4, index + 2) = DN_DX(iNode, 1);

            B(5, index + 0) = DN_DX(iNode, 2);
            B(5, index + 2) = DN_DX(iNode, 0);
        }
    }

    return B;

    KRATOS_CATCH("");
}

void StructuralMeshMovingElement::CheckElementMatrixDimension(MatrixType& rLeftHandSideMatrix,
                                                              VectorType& rRightHandSideVector)
{
    if (rLeftHandSideMatrix.size1() != mLocalSize)
        rLeftHandSideMatrix.resize(mLocalSize, mLocalSize, false);

    rLeftHandSideMatrix = ZeroMatrix(mLocalSize, mLocalSize);

    if (rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);
}

void StructuralMeshMovingElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                       VectorType& rRightHandSideVector,
                                                       ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    CheckElementMatrixDimension(rLeftHandSideMatrix, rRightHandSideVector);

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        double IntegrationWeight = integration_points[PointNumber].Weight();

        MatrixType B = CalculateBMatrix(dimension, PointNumber);

        MatrixType ConstitutiveMatrix = SetAndModifyConstitutiveLaw(dimension, PointNumber);
        // Compute LHS
        noalias(rLeftHandSideMatrix) +=
            prod(trans(B), IntegrationWeight * Matrix(prod(ConstitutiveMatrix, B)));

        // Compute RHS
        VectorType LastValues = ZeroVector(mLocalSize);
        this->GetDisplacementValues(LastValues, 0);
        noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, LastValues);
    }

    KRATOS_CATCH("");
}

void StructuralMeshMovingElement::EquationIdVector(EquationIdVectorType& rResult,
                                                   ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != mLocalSize)
        rResult.resize(mLocalSize, false);

    unsigned int pos = this->GetGeometry()[0].GetDofPosition(MESH_DISPLACEMENT_X);
    if (dimension == 2)
        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType Index = iNode * dimension;
            rResult[Index] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_X,pos).EquationId();
            rResult[Index + 1] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_Y,pos+1).EquationId();
        }
    else
        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType Index = iNode * dimension;
            rResult[Index] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_X,pos).EquationId();
            rResult[Index + 1] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_Y,pos+1).EquationId();
            rResult[Index + 2] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_Z,pos+2).EquationId();
        }
}

void StructuralMeshMovingElement::GetDofList(DofsVectorType& rElementalDofList,
                                             ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rElementalDofList.size() != mLocalSize)
        rElementalDofList.resize(mLocalSize);

    if (dimension == 2)
        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType Index = iNode * dimension;
            rElementalDofList[Index] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_X);
            rElementalDofList[Index + 1] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_Y);
        }
    else
        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType Index = iNode * dimension;
            rElementalDofList[Index] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_X);
            rElementalDofList[Index + 1] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_Y);
            rElementalDofList[Index + 2] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_Z);
        }
}

// Called in function "CalculateReactions" within the block builder and solver
void StructuralMeshMovingElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                         ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    MatrixType LeftHandSideMatrix = ZeroMatrix(mLocalSize, mLocalSize);

    if (rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);
    rRightHandSideVector.clear();

    // Get nodal solutions
    GeometryType& rGeom = this->GetGeometry();
    VectorType MeshDisplacement;
    MeshDisplacement.resize(mLocalSize, false);

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
    {
        MeshDisplacement[dimension * iNode + 0] =
            rGeom[iNode].GetSolutionStepValue(MESH_DISPLACEMENT_X);
        MeshDisplacement[dimension * iNode + 1] =
            rGeom[iNode].GetSolutionStepValue(MESH_DISPLACEMENT_Y);
        MeshDisplacement[dimension * iNode + 2] =
            rGeom[iNode].GetSolutionStepValue(MESH_DISPLACEMENT_Z);
    }

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        double IntegrationWeight = integration_points[PointNumber].Weight();

        MatrixType B = CalculateBMatrix(dimension, PointNumber);

        MatrixType ConstitutiveMatrix = SetAndModifyConstitutiveLaw(dimension, PointNumber);

        VectorType StressVector =
            prod(ConstitutiveMatrix, VectorType(prod(B, MeshDisplacement)));
        noalias(rRightHandSideVector) += IntegrationWeight * prod(trans(B), StressVector);
    }
}

// Only called in the CalculateLocalSystem Function above
void StructuralMeshMovingElement::CalculateAndAddRHS(VectorType& rRightHandSideVector)
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    if (rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);

    GeometryType& rGeom = this->GetGeometry();
    for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
    {
        // Note that we need to divide by the neighbours since the final RHS
        // needs
        // to be the desired value
        //(the addition below is done as many times as elements we have where
        // this node appears)

        int number_neighbours = rGeom[iNode].GetValue(NEIGHBOUR_ELEMENTS).size();

        rRightHandSideVector[dimension * iNode + 0] +=
            rGeom[iNode].GetSolutionStepValue(MESH_RHS_X) / number_neighbours;
        rRightHandSideVector[dimension * iNode + 1] +=
            rGeom[iNode].GetSolutionStepValue(MESH_RHS_Y) / number_neighbours;
        rRightHandSideVector[dimension * iNode + 2] +=
            rGeom[iNode].GetSolutionStepValue(MESH_RHS_Z) / number_neighbours;
    }
}
} // Namespace Kratos
