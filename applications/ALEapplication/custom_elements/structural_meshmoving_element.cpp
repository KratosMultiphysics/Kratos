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
    const GeometryType& r_geom = this->GetGeometry();

    return BaseType::Pointer(new StructuralMeshMovingElement(
        NewId, r_geom.Create(rThisNodes), pProperties));
}

Element::Pointer StructuralMeshMovingElement::Create(IndexType NewId,
                                                     GeometryType::Pointer pGeom,
                                                     PropertiesType::Pointer pProperties) const
{
    return BaseType::Pointer(new StructuralMeshMovingElement(NewId, pGeom, pProperties));
}

void StructuralMeshMovingElement::GetDisplacementValues(VectorType& rValues, const int Step)
{
    GeometryType& r_geom = this->GetGeometry();
    const SizeType num_nodes = r_geom.PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    if (dimension == 2)
    {
        SizeType index = 0;
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            rValues[index++] =
                r_geom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT_X, Step);
            rValues[index++] =
                r_geom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT_Y, Step);
        }
    }
    else if (dimension == 3)
    {
        SizeType index = 0;
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            rValues[index++] =
                r_geom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT_X, Step);
            rValues[index++] =
                r_geom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT_Y, Step);
            rValues[index++] =
                r_geom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT_Z, Step);
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
    mTotalDomainInitialSize = 0.0;

    const SizeType num_nodes = this->GetGeometry().PointsNumber();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    mLocalSize = num_nodes * dimension;

    for (unsigned int g = 0; g < integration_points.size(); ++g)
    {
        // getting informations for integration
        double weight = integration_points[g].Weight();

        // calculating and storing inverse of the jacobian and the parameters
        // needed
        MathUtils<double>::InvertMatrix(J0[g], mInvJ0[g], mDetJ0[g]);

        // calculating the total area
        mTotalDomainInitialSize += mDetJ0[g] * weight;
    }

    KRATOS_CATCH("");
}

StructuralMeshMovingElement::MatrixType StructuralMeshMovingElement::SetAndModifyConstitutiveLaw(
    const int& Dimension, const double& rPointNumber)
{
    KRATOS_TRY;

    // Stiffening of elements using Jacobian determinants and exponent between
    // 0.0 and 2.0
    const double J0 = 100; // Factor influences how far the displacement spreads
                           // into the fluid mesh
    const double xi = 1.5; // 1.5 Exponent influences stiffening of smaller
                           // elements; 0 = no stiffening
    const double quotient = J0 / mTotalDomainInitialSize;
    const double weight = mTotalDomainInitialSize * pow(quotient, xi);
    const double poisson_coefficient = 0.3;

    // The ratio between lamdbda and mu affects relative stiffening against
    // volume or shape change.
    double lambda = weight * poisson_coefficient /
                    ((1 + poisson_coefficient) * (1 - 2 * poisson_coefficient));
    double mu = weight / (2 * (1 - poisson_coefficient));

    MatrixType constitutive_matrix;

    // stress = lambda*tr(strain tensor)*I + 2*mu*(strain tensor).
    if (Dimension == 2)
    {
        constitutive_matrix = ZeroMatrix(3, 3);
        constitutive_matrix(0, 0) = lambda + 2 * mu;
        constitutive_matrix(1, 1) = constitutive_matrix(0, 0);
        constitutive_matrix(2, 2) = mu;
        constitutive_matrix(0, 1) = lambda;
        constitutive_matrix(1, 0) = lambda;
    }

    else if (Dimension == 3)
    {
        constitutive_matrix = ZeroMatrix(6, 6);
        constitutive_matrix(0, 0) = lambda + 2 * mu;
        constitutive_matrix(1, 1) = constitutive_matrix(0, 0);
        constitutive_matrix(2, 2) = constitutive_matrix(0, 0);
        constitutive_matrix(3, 3) = mu;
        constitutive_matrix(4, 4) = mu;
        constitutive_matrix(5, 5) = mu;
        constitutive_matrix(0, 1) = lambda;
        constitutive_matrix(1, 0) = lambda;
        constitutive_matrix(0, 2) = lambda;
        constitutive_matrix(2, 0) = lambda;
        constitutive_matrix(1, 2) = lambda;
        constitutive_matrix(2, 1) = lambda;
    }

    return constitutive_matrix;

    KRATOS_CATCH("");
}

StructuralMeshMovingElement::MatrixType StructuralMeshMovingElement::CalculateBMatrix(
    const int& Dimension, const double& rPointNumber)
{
    KRATOS_TRY;

    GeometryType::ShapeFunctionsGradientsType DN_De =
        this->GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    Matrix DN_DX = prod(DN_De[rPointNumber], mInvJ0[rPointNumber]);

    const SizeType num_nodes = this->GetGeometry().PointsNumber();

    MatrixType B;

    if (Dimension == 2)
    {
        B = ZeroMatrix(3, num_nodes * 2);

        SizeType index = 0;
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            B(0, index + 0) = DN_DX(i_node, 0);
            B(0, index + 1) = 0.0;
            B(1, index + 0) = 0.0;
            B(1, index + 1) = DN_DX(i_node, 1);
            B(2, index + 0) = DN_DX(i_node, 1);
            B(2, index + 1) = DN_DX(i_node, 0);
            index += 2;
        }
    }

    else if (Dimension == 3)
    {
        B = ZeroMatrix(6, num_nodes * 3);

        SizeType index = 0;
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            B(0, index + 0) = DN_DX(i_node, 0);
            B(1, index + 1) = DN_DX(i_node, 1);
            B(2, index + 2) = DN_DX(i_node, 2);
            B(3, index + 0) = DN_DX(i_node, 1);
            B(3, index + 1) = DN_DX(i_node, 0);
            B(4, index + 1) = DN_DX(i_node, 2);
            B(4, index + 2) = DN_DX(i_node, 1);
            B(5, index + 0) = DN_DX(i_node, 2);
            B(5, index + 2) = DN_DX(i_node, 0);
            index += 3;
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

    for (unsigned int g = 0; g < integration_points.size(); ++g)
    {
        double weight = integration_points[g].Weight();

        MatrixType B = CalculateBMatrix(dimension, g);

        MatrixType constitutive_matrix = SetAndModifyConstitutiveLaw(dimension, g);
        // Compute LHS
        noalias(rLeftHandSideMatrix) +=
            prod(trans(B), weight * Matrix(prod(constitutive_matrix, B)));

        // Compute RHS
        VectorType last_values;
        this->GetDisplacementValues(last_values, 0);
        noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, last_values);
    }

    KRATOS_CATCH("");
}

void StructuralMeshMovingElement::EquationIdVector(EquationIdVectorType& rResult,
                                                   ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& r_geom = this->GetGeometry();
    const SizeType num_nodes = r_geom.size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != mLocalSize)
        rResult.resize(mLocalSize, false);

    unsigned int pos = this->GetGeometry()[0].GetDofPosition(MESH_DISPLACEMENT_X);
    if (dimension == 2)
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            SizeType index = i_node * dimension;
            rResult[index] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_Y,pos+1).EquationId();
        }
    else
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            SizeType index = i_node * dimension;
            rResult[index] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_Y,pos+1).EquationId();
            rResult[index + 2] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_Z,pos+2).EquationId();
        }
}

void StructuralMeshMovingElement::GetDofList(DofsVectorType& rElementalDofList,
                                             ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& r_geom = this->GetGeometry();
    const SizeType num_nodes = r_geom.size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rElementalDofList.size() != mLocalSize)
        rElementalDofList.resize(mLocalSize);

    if (dimension == 2)
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            SizeType index = i_node * dimension;
            rElementalDofList[index] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_X);
            rElementalDofList[index + 1] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_Y);
        }
    else
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            SizeType index = i_node * dimension;
            rElementalDofList[index] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_X);
            rElementalDofList[index + 1] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_Y);
            rElementalDofList[index + 2] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_Z);
        }
}

// Called in function "CalculateReactions" within the block builder and solver
void StructuralMeshMovingElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                         ProcessInfo& rCurrentProcessInfo)
{
    MatrixType LHS;
    CalculateLocalSystem(LHS, rRightHandSideVector, rCurrentProcessInfo);
}

} // Namespace Kratos
