//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

// System includes

// External includes

// Project includes
#include "custom_elements/laplacian_meshmoving_element.h"
#include "includes/mesh_moving_variables.h"
#include "custom_utilities/move_mesh_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos {

LaplacianMeshMovingElement::LaplacianMeshMovingElement(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}
LaplacianMeshMovingElement::LaplacianMeshMovingElement(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

//******************************************************************************
//******************************************************************************
Element::Pointer
LaplacianMeshMovingElement::Create(IndexType NewId,
                                   NodesArrayType const &rThisNodes,
                                   PropertiesType::Pointer pProperties) const {
  return Kratos::make_shared<LaplacianMeshMovingElement>(
      NewId, GetGeometry().Create(rThisNodes), pProperties);
}

//******************************************************************************
//******************************************************************************
Element::Pointer
LaplacianMeshMovingElement::Create(IndexType NewId, GeometryType::Pointer pGeom,
                                   PropertiesType::Pointer pProperties) const {
  return Kratos::make_shared<LaplacianMeshMovingElement>(NewId, pGeom,
                                                         pProperties);
}

//******************************************************************************
//******************************************************************************
void LaplacianMeshMovingElement::CalculateDeltaPosition(
    VectorType &rIntermediateDisplacements, ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  const unsigned int component_index =
      rCurrentProcessInfo[LAPLACIAN_DIRECTION] - 1;
  const SizeType num_nodes = this->GetGeometry().PointsNumber();

  for (SizeType iNode = 0; iNode < num_nodes; ++iNode) {
    const VectorType &displacement =
        GetGeometry()[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT, 0) -
        GetGeometry()[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
    rIntermediateDisplacements[iNode] = displacement[component_index];
  }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void LaplacianMeshMovingElement::CheckElementMatrixDimension(
    MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector) {
  const SizeType num_nodes = this->GetGeometry().PointsNumber();

  if (rLeftHandSideMatrix.size1() != num_nodes)
    rLeftHandSideMatrix.resize(num_nodes, num_nodes, false);

  if (rRightHandSideVector.size() != num_nodes)
    rRightHandSideVector.resize(num_nodes, false);

  noalias(rLeftHandSideMatrix) = ZeroMatrix(num_nodes, num_nodes);
  noalias(rRightHandSideVector) = ZeroVector(num_nodes);
}

//******************************************************************************
//******************************************************************************
void LaplacianMeshMovingElement::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const IntegrationMethod this_integration_method =
      rgeom.GetDefaultIntegrationMethod();
  const GeometryType::IntegrationPointsArrayType &integration_points =
      GetGeometry().IntegrationPoints(this_integration_method);
  VectorType delta_displacement = ZeroVector(num_nodes);

  GeometryType::JacobiansType J0;
  GeometryType::JacobiansType invJ0;
  VectorType detJ0;

  MoveMeshUtilities::CheckJacobianDimension(invJ0, detJ0, rgeom);

  CheckElementMatrixDimension(rLeftHandSideMatrix, rRightHandSideVector);

  for (unsigned int point_number = 0; point_number < integration_points.size();
       ++point_number) {
    J0 = GetGeometry().Jacobian(J0, this_integration_method);
    MathUtils<double>::InvertMatrix(J0[point_number], invJ0[point_number],
                                    detJ0[point_number]);

    GeometryType::ShapeFunctionsGradientsType DN_De =
        this->GetGeometry().ShapeFunctionsLocalGradients(
            this_integration_method);

    // We do not multiply by DetJ0, since this additionally stabilizes the
    // simulation
    double integration_weight = integration_points[point_number].Weight();

    // Compute LHS
    Matrix DN_DX = prod(DN_De[point_number], invJ0[point_number]);
    noalias(rLeftHandSideMatrix) +=
        integration_weight * prod(DN_DX, trans(DN_DX));

    // Compute RHS
    CalculateDeltaPosition(delta_displacement, rCurrentProcessInfo);
    noalias(rRightHandSideVector) =
        -prod(rLeftHandSideMatrix, delta_displacement);
  }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void LaplacianMeshMovingElement::EquationIdVector(
    EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.size();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

  if (rResult.size() != num_nodes)
    rResult.resize(num_nodes, false);

  unsigned int pos = this->GetGeometry()[0].GetDofPosition(MESH_DISPLACEMENT_X);

  if (dimension == 2) {
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 1)
        rResult[i_node] =
            rgeom[i_node].GetDof(MESH_DISPLACEMENT_X, pos).EquationId();

      else if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 2)
        rResult[i_node] =
            rgeom[i_node].GetDof(MESH_DISPLACEMENT_Y, pos + 1).EquationId();
    }
  } else {
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 1)
        rResult[i_node] =
            rgeom[i_node].GetDof(MESH_DISPLACEMENT_X, pos).EquationId();
      else if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 2)
        rResult[i_node] =
            rgeom[i_node].GetDof(MESH_DISPLACEMENT_Y, pos + 1).EquationId();
      else if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 3)
        rResult[i_node] =
            rgeom[i_node].GetDof(MESH_DISPLACEMENT_Z, pos + 2).EquationId();
    }
  }
  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void LaplacianMeshMovingElement::GetDofList(DofsVectorType &rElementalDofList,
                                            ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.size();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

  if (rElementalDofList.size() != num_nodes)
    rElementalDofList.resize(num_nodes);

  if (dimension == 2)
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 1)
        rElementalDofList[i_node] = rgeom[i_node].pGetDof(MESH_DISPLACEMENT_X);
      else if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 2)
        rElementalDofList[i_node] = rgeom[i_node].pGetDof(MESH_DISPLACEMENT_Y);
    }
  else
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 1)
        rElementalDofList[i_node] = rgeom[i_node].pGetDof(MESH_DISPLACEMENT_X);
      if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 2)
        rElementalDofList[i_node] = rgeom[i_node].pGetDof(MESH_DISPLACEMENT_Y);
      if (rCurrentProcessInfo[LAPLACIAN_DIRECTION] == 3)
        rElementalDofList[i_node] = rgeom[i_node].pGetDof(MESH_DISPLACEMENT_Z);
    }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
// Called in function "CalculateReactions" within the component wise builder and
// solver
void LaplacianMeshMovingElement::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  MatrixType LHS;
  CalculateLocalSystem(LHS, rRightHandSideVector, rCurrentProcessInfo);

  KRATOS_CATCH("");
}

} // Namespace Kratos
