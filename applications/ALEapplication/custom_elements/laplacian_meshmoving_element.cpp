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

namespace Kratos
{
LaplacianMeshMovingElement::LaplacianMeshMovingElement(IndexType NewId,
                                                       GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

LaplacianMeshMovingElement::LaplacianMeshMovingElement(IndexType NewId,
                                                       GeometryType::Pointer pGeometry,
                                                       PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

Element::Pointer LaplacianMeshMovingElement::Create(IndexType NewId,
                                                    NodesArrayType const& rThisNodes,
                                                    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new LaplacianMeshMovingElement(
        NewId, GetGeometry().Create(rThisNodes), pProperties));
}

Element::Pointer LaplacianMeshMovingElement::Create(IndexType NewId,
                                                    GeometryType::Pointer pGeom,
                                                    PropertiesType::Pointer pProperties) const
{
    return BaseType::Pointer(new LaplacianMeshMovingElement(NewId, pGeom, pProperties));
}

void LaplacianMeshMovingElement::Initialize()
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    mJ0 = GetGeometry().Jacobian(mJ0, mThisIntegrationMethod);
    mInvJ0.resize(integration_points.size());
    mDetJ0.resize(integration_points.size(), false);
    mTotalDomainInitialSize = 0.00;

    const SizeType num_nodes = this->GetGeometry().PointsNumber();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    mLocalSize = num_nodes * dimension;

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        // getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();

        // calculating and storing inverse of the jacobian and the parameters
        // needed
        MathUtils<double>::InvertMatrix(
            mJ0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber]);

        // calculating the total area
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }

    KRATOS_CATCH("");
}


void LaplacianMeshMovingElement::GetDisplacementValues(VectorType& rValues, const int Step)
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

LaplacianMeshMovingElement::MatrixType LaplacianMeshMovingElement::CalculateDerivatives(
    const int& rdimension, const double& rPointNumber)
{
    KRATOS_TRY;

    GeometryType::ShapeFunctionsGradientsType DN_De =
        this->GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    Matrix DN_DX = prod(DN_De[rPointNumber], mInvJ0[rPointNumber]);

    return DN_DX;

    KRATOS_CATCH("");
}

void LaplacianMeshMovingElement::CalculateDeltaPosition(VectorType& rIntermediateDisplacements,
                                                        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const unsigned int ComponentIndex = rCurrentProcessInfo[FRACTIONAL_STEP] - 1;

    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
    {
        const VectorType& displacement =
            GetGeometry()[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT, 0) -
            GetGeometry()[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
        rIntermediateDisplacements[iNode] = displacement[ComponentIndex];
    }

    KRATOS_CATCH("");
}

void LaplacianMeshMovingElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                      VectorType& rRightHandSideVector,
                                                      ProcessInfo& rCurrentProcessInfo)
{
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

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        double IntegrationWeight = integration_points[PointNumber].Weight() * mDetJ0[PointNumber];

        // Compute LHS
        MatrixType DN_DX = CalculateDerivatives(dimension, PointNumber);
        noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(DN_DX, trans(DN_DX));

    }

    Vector temp(NumNodes);
    // Compute RHS
    double LocalSize = NumNodes * dimension;
    VectorType last_values = ZeroVector(LocalSize);

    unsigned int ComponentIndex = rCurrentProcessInfo[FRACTIONAL_STEP] - 1;

  const array_1d<double, 3>& disp0 = GetGeometry()[0].FastGetSolutionStepValue(
      MESH_DISPLACEMENT, 0)
      - GetGeometry()[0].FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
  const array_1d<double, 3>& disp1 = GetGeometry()[1].FastGetSolutionStepValue(
      MESH_DISPLACEMENT, 0)
      - GetGeometry()[1].FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
  const array_1d<double, 3>& disp2 = GetGeometry()[2].FastGetSolutionStepValue(
      MESH_DISPLACEMENT, 0)
      - GetGeometry()[2].FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);



  array_1d<double, 3> temp_vec_np;
  //VectorType temp_vec_np;

  temp_vec_np[0] = disp0[ComponentIndex];
  temp_vec_np[1] = disp1[ComponentIndex];
  temp_vec_np[2] = disp2[ComponentIndex];


  

  KRATOS_WATCH(ComponentIndex);


    this->GetDisplacementValues(last_values, 1);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, temp_vec_np);
    
    KRATOS_CATCH("");
}

void LaplacianMeshMovingElement::EquationIdVector(EquationIdVectorType& rResult,
                                                   ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& r_geom = this->GetGeometry();
    const SizeType num_nodes = r_geom.size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != num_nodes)
        rResult.resize(num_nodes, false);

    unsigned int pos = this->GetGeometry()[0].GetDofPosition(MESH_DISPLACEMENT_X);

    if (dimension == 2)
        {
            for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
            {
                //SizeType index = i_node * dimension;

                if(rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
                rResult[i_node] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_X,pos).EquationId();

                else if(rCurrentProcessInfo[FRACTIONAL_STEP] == 2)
                rResult[i_node] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_Y,pos+1).EquationId();
            }
        }
    else
        {
            for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
            {
                //SizeType index = i_node * dimension;

                if(rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
                rResult[i_node] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_X,pos).EquationId();
                else if(rCurrentProcessInfo[FRACTIONAL_STEP] == 2)
                rResult[i_node] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_Y,pos+1).EquationId();
                else if(rCurrentProcessInfo[FRACTIONAL_STEP] == 3)
                rResult[i_node] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_Z,pos+2).EquationId();
            }
        }
}

void LaplacianMeshMovingElement::GetDofList(DofsVectorType& rElementalDofList,
                                             ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& r_geom = this->GetGeometry();
    const SizeType num_nodes = r_geom.size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rElementalDofList.size() != num_nodes)
        rElementalDofList.resize(num_nodes);

    if (dimension == 2)
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            //SizeType index = i_node * dimension;

            if(rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
            rElementalDofList[i_node] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_X);
            else if (rCurrentProcessInfo[FRACTIONAL_STEP] == 2)
            rElementalDofList[i_node] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_Y);
        }
    else
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            //SizeType index = i_node * dimension;

            if(rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
            rElementalDofList[i_node] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_X);
            if(rCurrentProcessInfo[FRACTIONAL_STEP] == 2)
            rElementalDofList[i_node] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_Y);
            if(rCurrentProcessInfo[FRACTIONAL_STEP] == 3)
            rElementalDofList[i_node] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_Z);
        }
}

} // Namespace Kratos
