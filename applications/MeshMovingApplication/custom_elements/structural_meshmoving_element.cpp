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
//#include "includes/define.h"
#include "custom_elements/structural_meshmoving_element.h"
#include "includes/mesh_moving_variables.h"
#include "custom_utilities/move_mesh_utilities.h"

namespace Kratos {
StructuralMeshMovingElement::StructuralMeshMovingElement(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

StructuralMeshMovingElement::StructuralMeshMovingElement(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

//******************************************************************************
//******************************************************************************
Element::Pointer
StructuralMeshMovingElement::Create(IndexType NewId,
                                    NodesArrayType const &rThisNodes,
                                    PropertiesType::Pointer pProperties) const {
  const GeometryType &rgeom = this->GetGeometry();

  return Kratos::make_shared<StructuralMeshMovingElement>(
      NewId, rgeom.Create(rThisNodes), pProperties);
}

//******************************************************************************
//******************************************************************************
Element::Pointer
StructuralMeshMovingElement::Create(IndexType NewId,
                                    GeometryType::Pointer pGeom,
                                    PropertiesType::Pointer pProperties) const {
  return Kratos::make_shared<StructuralMeshMovingElement>(NewId, pGeom,
                                                          pProperties);
}


//******************************************************************************
//******************************************************************************
void StructuralMeshMovingElement::GetValuesVector(VectorType &rValues,
                                                  int Step) {
  GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rValues.size() != local_size)
    rValues.resize(local_size, false);

  if (dimension == 2) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT_X, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT_Y, Step);
    }
  } else if (dimension == 3) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT_X, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT_Y, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT_Z, Step);
    }
  }
}

//******************************************************************************
//******************************************************************************
StructuralMeshMovingElement::MatrixType
StructuralMeshMovingElement::SetAndModifyConstitutiveLaw(
    const int &Dimension, const double &rPointNumber) {
  KRATOS_TRY;

  GeometryType::JacobiansType J0;
  GeometryType::JacobiansType invJ0;
  VectorType detJ0;
  GeometryType &rgeom = this->GetGeometry();
  const IntegrationMethod this_integration_method =
      rgeom.GetDefaultIntegrationMethod();

  MoveMeshUtilities::CheckJacobianDimension(invJ0, detJ0, rgeom);

  J0 = GetGeometry().Jacobian(J0, this_integration_method);

  MathUtils<double>::InvertMatrix(J0[rPointNumber], invJ0[rPointNumber],
                                  detJ0[rPointNumber]);

  // Stiffening of elements using Jacobian determinants and exponent between
  // 0.0 and 2.0
  const double factor =
      100;               // Factor influences how far the displacement spreads
                         // into the fluid mesh
  const double xi = 1.5; // 1.5 Exponent influences stiffening of smaller
                         // elements; 0 = no stiffening
  const double quotient = factor / detJ0[rPointNumber];
  const double weighting_factor = detJ0[rPointNumber] * std::pow(quotient, xi);
  const double poisson_coefficient = 0.3;

  // The ratio between lambda and mu affects relative stiffening against
  // volume or shape change.
  const double lambda =
      weighting_factor * poisson_coefficient /
      ((1 + poisson_coefficient) * (1 - 2 * poisson_coefficient));
  const double mu = weighting_factor / (2 * (1 + poisson_coefficient));

  MatrixType constitutive_matrix;

  // stress = lambda*tr(strain tensor)*I + 2*mu*(strain tensor).
  if (Dimension == 2) {
    constitutive_matrix = ZeroMatrix(3, 3);
    constitutive_matrix(0, 0) = lambda + 2 * mu;
    constitutive_matrix(1, 1) = constitutive_matrix(0, 0);
    constitutive_matrix(2, 2) = mu;
    constitutive_matrix(0, 1) = lambda;
    constitutive_matrix(1, 0) = lambda;
  }

  else if (Dimension == 3) {
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

//******************************************************************************
//******************************************************************************
StructuralMeshMovingElement::MatrixType
StructuralMeshMovingElement::CalculateBMatrix(const int &Dimension,
                                              const double &rPointNumber) {
  KRATOS_TRY;

  GeometryType &rgeom = this->GetGeometry();
  const IntegrationMethod this_integration_method =
      rgeom.GetDefaultIntegrationMethod();
  GeometryType::ShapeFunctionsGradientsType DN_De =
      rgeom.ShapeFunctionsLocalGradients(this_integration_method);
  GeometryType::JacobiansType J0;
  GeometryType::JacobiansType invJ0;
  VectorType detJ0;

  MoveMeshUtilities::CheckJacobianDimension(invJ0, detJ0, rgeom);

  J0 = GetGeometry().Jacobian(J0, this_integration_method);
  MathUtils<double>::InvertMatrix(J0[rPointNumber], invJ0[rPointNumber],
                                  detJ0[rPointNumber]);

  Matrix DN_DX = prod(DN_De[rPointNumber], invJ0[rPointNumber]);

  const SizeType num_nodes = rgeom.PointsNumber();

  MatrixType B;

  if (Dimension == 2) {
    B = ZeroMatrix(3, num_nodes * 2);

    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      B(0, index + 0) = DN_DX(i_node, 0);
      B(0, index + 1) = 0.0;
      B(1, index + 0) = 0.0;
      B(1, index + 1) = DN_DX(i_node, 1);
      B(2, index + 0) = DN_DX(i_node, 1);
      B(2, index + 1) = DN_DX(i_node, 0);
      index += 2;
    }
  }

  else if (Dimension == 3) {
    B = ZeroMatrix(6, num_nodes * 3);

    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
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

//******************************************************************************
//******************************************************************************
void StructuralMeshMovingElement::CheckElementMatrixDimension(
    MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector) {
  GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const unsigned int dimension = rgeom.WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rLeftHandSideMatrix.size1() != local_size)
    rLeftHandSideMatrix.resize(local_size, local_size, false);

  noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);

  if (rRightHandSideVector.size() != local_size)
    rRightHandSideVector.resize(local_size, false);
}

//******************************************************************************
//******************************************************************************
void StructuralMeshMovingElement::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  GeometryType &rgeom = this->GetGeometry();
  const IntegrationMethod this_integration_method =
      rgeom.GetDefaultIntegrationMethod();
  const unsigned int dimension = rgeom.WorkingSpaceDimension();

  const GeometryType::IntegrationPointsArrayType &integration_points =
      GetGeometry().IntegrationPoints(this_integration_method);

  CheckElementMatrixDimension(rLeftHandSideMatrix, rRightHandSideVector);

  for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number) {
    double weight = integration_points[point_number].Weight();

    MatrixType B = CalculateBMatrix(dimension, point_number);

    MatrixType constitutive_matrix = SetAndModifyConstitutiveLaw(dimension, point_number);
    // Compute LHS
    noalias(rLeftHandSideMatrix) +=
        prod(trans(B), weight * Matrix(prod(constitutive_matrix, B)));

    // Compute RHS
    VectorType last_values;
    this->GetValuesVector(last_values, 0);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, last_values);
  }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void StructuralMeshMovingElement::EquationIdVector(
    EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.size();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rResult.size() != local_size)
    rResult.resize(local_size, false);

  unsigned int pos = rgeom[0].GetDofPosition(MESH_DISPLACEMENT_X);
  if (dimension == 2)
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      SizeType index = i_node * dimension;
      rResult[index] =
          rgeom[i_node].GetDof(MESH_DISPLACEMENT_X, pos).EquationId();
      rResult[index + 1] =
          rgeom[i_node].GetDof(MESH_DISPLACEMENT_Y, pos + 1).EquationId();
    }
  else
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      SizeType index = i_node * dimension;
      rResult[index] =
          rgeom[i_node].GetDof(MESH_DISPLACEMENT_X, pos).EquationId();
      rResult[index + 1] =
          rgeom[i_node].GetDof(MESH_DISPLACEMENT_Y, pos + 1).EquationId();
      rResult[index + 2] =
          rgeom[i_node].GetDof(MESH_DISPLACEMENT_Z, pos + 2).EquationId();
    }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void StructuralMeshMovingElement::GetDofList(DofsVectorType &rElementalDofList,
                                             ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.size();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rElementalDofList.size() != local_size)
    rElementalDofList.resize(local_size);

  if (dimension == 2)
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      SizeType index = i_node * dimension;
      rElementalDofList[index] = rgeom[i_node].pGetDof(MESH_DISPLACEMENT_X);
      rElementalDofList[index + 1] = rgeom[i_node].pGetDof(MESH_DISPLACEMENT_Y);
    }
  else
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      SizeType index = i_node * dimension;
      rElementalDofList[index] = rgeom[i_node].pGetDof(MESH_DISPLACEMENT_X);
      rElementalDofList[index + 1] = rgeom[i_node].pGetDof(MESH_DISPLACEMENT_Y);
      rElementalDofList[index + 2] = rgeom[i_node].pGetDof(MESH_DISPLACEMENT_Z);
    }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
// Called in function "CalculateReactions" within the block builder and solver
void StructuralMeshMovingElement::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  MatrixType LHS;
  CalculateLocalSystem(LHS, rRightHandSideVector, rCurrentProcessInfo);

  KRATOS_CATCH("");
}

} // Namespace Kratos
