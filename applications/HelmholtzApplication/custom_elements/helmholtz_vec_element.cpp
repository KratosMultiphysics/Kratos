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
//  Main authors:    Reza Najian Asl
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "helmholtz_application_variables.h"
#include "custom_elements/helmholtz_vec_element.h"

namespace Kratos {
HelmholtzVecElement::HelmholtzVecElement(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

HelmholtzVecElement::HelmholtzVecElement(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

//******************************************************************************
//******************************************************************************
Element::Pointer
HelmholtzVecElement::Create(IndexType NewId,
                                    NodesArrayType const &rThisNodes,
                                    PropertiesType::Pointer pProperties) const {
  const GeometryType &rgeom = this->GetGeometry();

  return Kratos::make_intrusive<HelmholtzVecElement>(
      NewId, rgeom.Create(rThisNodes), pProperties);
}

//******************************************************************************
//******************************************************************************
Element::Pointer
HelmholtzVecElement::Create(IndexType NewId,
                                    GeometryType::Pointer pGeom,
                                    PropertiesType::Pointer pProperties) const {
  return Kratos::make_intrusive<HelmholtzVecElement>(NewId, pGeom,
                                                          pProperties);
}


//******************************************************************************
//******************************************************************************
void HelmholtzVecElement::GetValuesVector(VectorType &rValues,
                                                  int Step) const {
  const GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rValues.size() != local_size)
    rValues.resize(local_size, false);

  if (dimension == 2) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_X, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Y, Step);
    }
  } else if (dimension == 3) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_X, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Y, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Z, Step);
    }
  }
}

//******************************************************************************
//******************************************************************************
HelmholtzVecElement::MatrixType
HelmholtzVecElement::SetAndModifyConstitutiveLaw(
    const int Dimension, const int PointNumber) const {
  KRATOS_TRY;

  GeometryType::JacobiansType J0;
  GeometryType::JacobiansType invJ0;
  VectorType detJ0;
  const GeometryType &rgeom = this->GetGeometry();
  const IntegrationMethod this_integration_method =
      rgeom.GetDefaultIntegrationMethod();


  const auto& r_integration_points = rgeom.IntegrationPoints(this_integration_method);
  if (invJ0.size() != r_integration_points.size()) {
    invJ0.resize(r_integration_points.size());
  }
  if (detJ0.size() != r_integration_points.size()) {
    detJ0.resize(r_integration_points.size());
  }

  J0 = GetGeometry().Jacobian(J0, this_integration_method);

  MathUtils<double>::InvertMatrix(J0[PointNumber], invJ0[PointNumber],
                                  detJ0[PointNumber]);

  // Stiffening of elements using Jacobian determinants and exponent between
  // 0.0 and 2.0
  const double factor =
      100;               // Factor influences how far the displacement spreads
                         // into the fluid mesh
  const double xi = 1.5; // 1.5 Exponent influences stiffening of smaller
                         // elements; 0 = no stiffening
  const double quotient = factor / detJ0[PointNumber];
  const double weighting_factor = detJ0[PointNumber] * std::pow(quotient, xi);
  const double poisson_coefficient = this->pGetProperties()->Has(HELMHOLTZ_POISSON_RATIO)
    ? this->pGetProperties()->GetValue(HELMHOLTZ_POISSON_RATIO) : 0.3;

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
HelmholtzVecElement::MatrixType
HelmholtzVecElement::CalculateBMatrix(const int Dimension,
                                              const int PointNumber) const {
  KRATOS_TRY;

  const GeometryType &rgeom = this->GetGeometry();
  const IntegrationMethod this_integration_method =
      rgeom.GetDefaultIntegrationMethod();
  GeometryType::ShapeFunctionsGradientsType DN_De =
      rgeom.ShapeFunctionsLocalGradients(this_integration_method);
  GeometryType::JacobiansType J0;
  GeometryType::JacobiansType invJ0;
  VectorType detJ0;

  const auto& r_integration_points = rgeom.IntegrationPoints(this_integration_method);
  if (invJ0.size() != r_integration_points.size()) {
    invJ0.resize(r_integration_points.size());
  }
  if (detJ0.size() != r_integration_points.size()) {
    detJ0.resize(r_integration_points.size());
  }

  J0 = GetGeometry().Jacobian(J0, this_integration_method);
  MathUtils<double>::InvertMatrix(J0[PointNumber], invJ0[PointNumber],
                                  detJ0[PointNumber]);

  Matrix DN_DX = prod(DN_De[PointNumber], invJ0[PointNumber]);

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
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVecElement::CalculateMMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix( mat_size, mat_size );

    // CONSISTENT MASS

    Element::GeometryType::JacobiansType J0;
    r_geom.Jacobian(J0,r_geom.GetDefaultIntegrationMethod());
    Matrix InvJ0(dimension,dimension);
    double detJ0;


    const IntegrationMethod integration_method = r_geom.GetDefaultIntegrationMethod();
    const auto& integration_points = r_geom.IntegrationPoints(integration_method);
    const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        //calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[point_number],InvJ0,detJ0);;
        const double integration_weight = integration_points[point_number].Weight() * detJ0;
        
        const Vector& rN = row(Ncontainer,point_number);

        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            const SizeType index_i = i * dimension;

            for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                const SizeType index_j = j * dimension;
                const double NiNj_weight = rN[i] * rN[j] * integration_weight;

                for ( IndexType k = 0; k < dimension; ++k )
                    rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
            }
        }
    }

    KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void HelmholtzVecElement::CheckElementMatrixDimension(
    MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector) const {
  const GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const unsigned int dimension = rgeom.WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rLeftHandSideMatrix.size1() != local_size)
    rLeftHandSideMatrix.resize(local_size, local_size, false);

  noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);

  if (rRightHandSideVector.size() != local_size)
    rRightHandSideVector.resize(local_size, false);

  noalias(rRightHandSideVector) = ZeroVector(local_size);
}

//******************************************************************************
//******************************************************************************
void HelmholtzVecElement::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  GeometryType &rgeom = this->GetGeometry();
  const IntegrationMethod this_integration_method =
      rgeom.GetDefaultIntegrationMethod();
  const unsigned int dimension = rgeom.WorkingSpaceDimension();

  const GeometryType::IntegrationPointsArrayType &integration_points =
      GetGeometry().IntegrationPoints(this_integration_method);

  CheckElementMatrixDimension(rLeftHandSideMatrix, rRightHandSideVector);


    Element::GeometryType::JacobiansType J0;
    rgeom.Jacobian(J0,rgeom.GetDefaultIntegrationMethod());
    Matrix InvJ0(dimension,dimension);
    double detJ0;

  for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number) {
     MathUtils<double>::InvertMatrix(J0[point_number],InvJ0,detJ0);;
    double weight = integration_points[point_number].Weight() * detJ0;

    MatrixType B = CalculateBMatrix(dimension, point_number);

    MatrixType constitutive_matrix = SetAndModifyConstitutiveLaw(dimension, point_number);
    // Compute LHS
    noalias(rLeftHandSideMatrix) +=
        -400 * prod(trans(B), weight * Matrix(prod(constitutive_matrix, B)));

    // Compute RHS
    VectorType last_values;
    this->GetValuesVector(last_values, 0);
    // noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, last_values);
  }

  MatrixType M;
  CalculateMMatrix(M,rCurrentProcessInfo);

  const unsigned int number_of_points = rgeom.size();
  Vector nodal_vals(number_of_points*3);
  for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
  {
      const VectorType &source = rgeom[node_element].FastGetSolutionStepValue(HELMHOLTZ_SOURCE);
      nodal_vals[3 * node_element + 0] = source[0];
      nodal_vals[3 * node_element + 1] = source[1];
      nodal_vals[3 * node_element + 2] = source[2];
  }  

  noalias(rLeftHandSideMatrix) += M;

  noalias(rRightHandSideVector) += prod(M,nodal_vals);

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void HelmholtzVecElement::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const {
  KRATOS_TRY;

  const GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.size();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rResult.size() != local_size)
    rResult.resize(local_size, false);

  unsigned int pos = rgeom[0].GetDofPosition(HELMHOLTZ_VARS_X);
  if (dimension == 2)
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      SizeType index = i_node * dimension;
      rResult[index] =
          rgeom[i_node].GetDof(HELMHOLTZ_VARS_X, pos).EquationId();
      rResult[index + 1] =
          rgeom[i_node].GetDof(HELMHOLTZ_VARS_Y, pos + 1).EquationId();
    }
  else
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      SizeType index = i_node * dimension;
      rResult[index] =
          rgeom[i_node].GetDof(HELMHOLTZ_VARS_X, pos).EquationId();
      rResult[index + 1] =
          rgeom[i_node].GetDof(HELMHOLTZ_VARS_Y, pos + 1).EquationId();
      rResult[index + 2] =
          rgeom[i_node].GetDof(HELMHOLTZ_VARS_Z, pos + 2).EquationId();
    }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void HelmholtzVecElement::GetDofList(DofsVectorType &rElementalDofList,
                                             const ProcessInfo &rCurrentProcessInfo) const {
  KRATOS_TRY;

  const GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.size();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rElementalDofList.size() != local_size)
    rElementalDofList.resize(local_size);

  if (dimension == 2)
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      SizeType index = i_node * dimension;
      rElementalDofList[index] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_X);
      rElementalDofList[index + 1] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_Y);
    }
  else
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      SizeType index = i_node * dimension;
      rElementalDofList[index] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_X);
      rElementalDofList[index + 1] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_Y);
      rElementalDofList[index + 2] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_Z);
    }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
// Called in function "CalculateReactions" within the block builder and solver
void HelmholtzVecElement::CalculateRightHandSide(
    VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  MatrixType LHS;
  CalculateLocalSystem(LHS, rRightHandSideVector, rCurrentProcessInfo);

  KRATOS_CATCH("");
}

int HelmholtzVecElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    Element::Check(rCurrentProcessInfo);

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( const auto& r_node : GetGeometry() ) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HELMHOLTZ_VARS,r_node)

        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_Z, r_node)
    }

    return 0;

    KRATOS_CATCH( "" );
}

} // Namespace Kratos
