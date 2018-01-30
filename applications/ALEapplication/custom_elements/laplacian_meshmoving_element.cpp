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
 *  Date:                $Date: February 2018 $
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
}

LaplacianMeshMovingElement::LaplacianMeshMovingElement(IndexType NewId,
                                                       GeometryType::Pointer pGeometry,
                                                       PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
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


void LaplacianMeshMovingElement::CalculateInitialJacobianValues(MatrixType& rJ0, MatrixType& rInvJ0, double& rDetJ0, const double& rPointNumber)
{
    KRATOS_TRY;

    const IntegrationMethod ThisIntegrationMethod = this->GetGeometry().GetDefaultIntegrationMethod();
    const MatrixType& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[rPointNumber];

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        const array_1d<double, 3>& coords = GetGeometry()[i].GetInitialPosition(); //NOTE: here we refer to the original, undeformed position!!
        for(unsigned int k = 0; k < GetGeometry().WorkingSpaceDimension(); k++)
        {
            for(unsigned int m = 0; m < GetGeometry().LocalSpaceDimension(); m++)
            {
                rJ0(k, m) += coords[k]*DN_De(i, m);
            }
        }
    }

    MathUtils<double>::InvertMatrix(rJ0, rInvJ0, rDetJ0);

    KRATOS_CATCH("");
}


LaplacianMeshMovingElement::MatrixType LaplacianMeshMovingElement::CalculateDerivatives(
    const int& rdimension, const double& rPointNumber)
{
    KRATOS_TRY;

    const IntegrationMethod ThisIntegrationMethod = this->GetGeometry().GetDefaultIntegrationMethod();
    uint dimension = GetGeometry().WorkingSpaceDimension();
    uint number_of_nodes = GetGeometry().size();
    MatrixType DN_DX( number_of_nodes, dimension);
    MatrixType rJ0(dimension, dimension);
    MatrixType rinvJ0(dimension, dimension);

    rJ0.clear();
        
    double rdetJ0;    
    const MatrixType& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[rPointNumber];

    CalculateInitialJacobianValues(rJ0, rinvJ0, rdetJ0, rPointNumber);

    noalias(DN_DX) = prod(DN_De, rinvJ0);

    return DN_DX;

    KRATOS_CATCH("");
}

void LaplacianMeshMovingElement::CalculateDeltaPosition(VectorType& rIntermediateDisplacements,
                                                        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const unsigned int ComponentIndex = rCurrentProcessInfo[LAPLACIAN_DIRECTION] - 1;
    const SizeType num_nodes = this->GetGeometry().PointsNumber();

    for (SizeType iNode = 0; iNode < num_nodes; ++iNode)
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

    const SizeType num_nodes = this->GetGeometry().PointsNumber();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    const IntegrationMethod ThisIntegrationMethod = this->GetGeometry().GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(ThisIntegrationMethod);

    if (rLeftHandSideMatrix.size1() != num_nodes)
        rLeftHandSideMatrix.resize(num_nodes, num_nodes, false);

    if (rRightHandSideVector.size() != num_nodes)
        rRightHandSideVector.resize(num_nodes, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(num_nodes, num_nodes);
    noalias(rRightHandSideVector) = ZeroVector(num_nodes);

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        // We do not multiply by DetJ0 since this stabilizes the simulation
        double IntegrationWeight = integration_points[PointNumber].Weight(); 

        // Compute LHS
        Matrix DN_DX = CalculateDerivatives(dimension, PointNumber);
        noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(DN_DX, trans(DN_DX));

    }

    // Compute RHS
    VectorType delta_displacement = ZeroVector(3);
    CalculateDeltaPosition(delta_displacement, rCurrentProcessInfo);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, delta_displacement);
    
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
                if(rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 1)
                rResult[i_node] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_X,pos).EquationId();

                else if(rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 2)
                rResult[i_node] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_Y,pos+1).EquationId();
            }
        }
    else
        {
            for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
            {
                if(rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 1)
                rResult[i_node] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_X,pos).EquationId();
                else if(rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 2)
                rResult[i_node] = r_geom[i_node].GetDof(MESH_DISPLACEMENT_Y,pos+1).EquationId();
                else if(rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 3)
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
            if(rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 1)
            rElementalDofList[i_node] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_X);
            else if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 2)
            rElementalDofList[i_node] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_Y);
        }
    else
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
        {
            if(rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 1)
            rElementalDofList[i_node] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_X);
            if(rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 2)
            rElementalDofList[i_node] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_Y);
            if(rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 3)
            rElementalDofList[i_node] = r_geom[i_node].pGetDof(MESH_DISPLACEMENT_Z);
        }
}

} // Namespace Kratos
